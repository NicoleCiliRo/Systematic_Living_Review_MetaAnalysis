"""
meta_niri.py

Meta-análisis de exactitud diagnóstica (NIRI para caries interproximal)
usando Python + rpy2 + R(mada::reitsma).

Requiere:
- R instalado
- Paquetes Python: rpy2, pandas, matplotlib
- Paquete R: mada (se instala automáticamente si no está)

Estructura mínima del CSV:
study_id,setting,sensitivity,specificity,auc,n_diseased,n_non_diseased,...

"""

import os
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector

pandas2ri.activate()

# -------------------------------------------------------------------
# 1. Configuración básica
# -------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
DATA_PATH = BASE_DIR / "data" / "niri_meta.csv"
OUTPUT_DIR = BASE_DIR / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)


# -------------------------------------------------------------------
# 2. Setup de R y paquete mada
# -------------------------------------------------------------------

def setup_r_mada():
    """Asegura que el paquete 'mada' está instalado en R y lo importa."""
    utils = importr('utils')

    # Elegir un mirror de CRAN (solo primera vez; luego puedes comentar)
    try:
        utils.chooseCRANmirror(ind=1)
    except Exception:
        # Si ya hay mirror elegido o hay algún problema, seguimos.
        pass

    installed = utils.installed_packages()
    installed_names = list(installed.rx(True, 'Package'))

    if 'mada' not in installed_names:
        utils.install_packages(StrVector(['mada']))

    mada = importr('mada')
    base = importr('base')
    return mada, base


# -------------------------------------------------------------------
# 3. Carga de datos y reconstrucción de TP/FN/FP/TN
# -------------------------------------------------------------------

def load_and_prepare_data(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)

    # Renombrar por si las columnas vienen con otros nombres
    rename_map = {
        'Study_ID': 'study_id',
        'Study': 'study_id',
        'Setting': 'setting',
        'Sens': 'sensitivity',
        'Spec': 'specificity',
        'AUC': 'auc',
        'N_pos': 'n_diseased',
        'N_neg': 'n_non_diseased'
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    required = ['study_id', 'setting', 'sensitivity', 'specificity',
                'n_diseased', 'n_non_diseased']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Faltan columnas necesarias en el CSV: {missing}")

    # Construir 2x2 aproximado a partir de Se, Sp y tamaños.
    # Redondeamos al entero más cercano.
    df['TP'] = (df['sensitivity'] * df['n_diseased']).round().astype(int)
    df['FN'] = df['n_diseased'] - df['TP']
    df['TN'] = (df['specificity'] * df['n_non_diseased']).round().astype(int)
    df['FP'] = df['n_non_diseased'] - df['TN']

    # Cálculo de Se/Sp "recomputadas" y sus errores estándar (para forest plots)
    df['sens_recalc'] = df['TP'] / (df['TP'] + df['FN'])
    df['spec_recalc'] = df['TN'] / (df['TN'] + df['FP'])

    df['se_sens'] = (df['sens_recalc'] * (1 - df['sens_recalc']) /
                     (df['TP'] + df['FN'])).pow(0.5)
    df['se_spec'] = (df['spec_recalc'] * (1 - df['spec_recalc']) /
                     (df['TN'] + df['FP'])).pow(0.5)

    return df


# -------------------------------------------------------------------
# 4. Función auxiliar: ejecutar Reitsma (bivariado) vía mada
# -------------------------------------------------------------------

def run_reitsma(mada, r_df, moderators=None):
    """
    Ejecuta el modelo Reitsma en R (mada::reitsma).

    Parameters
    ----------
    mada : R package object
    r_df : R data.frame
    moderators : list[str] or None
        Lista de nombres de covariables para meta-regresión, p.ej. ['setting', 'scanner'].

    Returns
    -------
    fit : R object
    summary_str : str
    """
    r_summary = ro.r['summary']

    if moderators:
        formula_str = "~ " + " + ".join(moderators)
        formula = ro.r(formula_str)
        fit = mada.reitsma(
            data=r_df,
            TP="TP",
            FN="FN",
            FP="FP",
            TN="TN",
            formula=formula
        )
    else:
        fit = mada.reitsma(
            data=r_df,
            TP="TP",
            FN="FN",
            FP="FP",
            TN="TN"
        )

    # Capturamos el summary como texto para guardarlo
    capture_output = ro.r('capture.output')
    summary_out = capture_output(r_summary(fit))
    summary_str = "\n".join(list(summary_out))

    return fit, summary_str


# -------------------------------------------------------------------
# 5. SROC desde mada -> pandas -> plot en Python
# -------------------------------------------------------------------

def get_sroc_df(mada, fit) -> pd.DataFrame:
    r_sroc = mada.sroc(fit)
    sroc_df = pandas2ri.rpy2py(r_sroc)
    # mada::sroc devuelve típicamente fpr y sens
    if 'fpr' not in sroc_df.columns:
        # Por si el nombre cambia, intentamos localizar algo similar
        raise RuntimeError("No encuentro columna 'fpr' en la salida de sroc(). Revisa la versión de mada.")
    sroc_df['specificity'] = 1 - sroc_df['fpr']
    return sroc_df


def plot_sroc(sroc_df: pd.DataFrame, title: str, filename: Path):
    plt.figure(figsize=(6, 6), dpi=300)
    plt.plot(sroc_df['fpr'], sroc_df['sens'], label='SROC')
    plt.plot([0, 1], [0, 1], linestyle='--', label='No discrimination')
    plt.xlabel('False Positive Rate (1 - Specificity)')
    plt.ylabel('Sensitivity')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


# -------------------------------------------------------------------
# 6. Forest plots (Sensibilidad, Especificidad, AUC)
# -------------------------------------------------------------------

def forest_plot_metric(df: pd.DataFrame,
                       metric: str,
                       se_metric: str | None,
                       title: str,
                       xlabel: str,
                       filename: Path):
    """
    Forest plot simple para una métrica (sens, spec, auc).

    Si se_metric es None, no pinta barras de error (solo puntos).
    """
    df = df.copy().reset_index(drop=True)
    y = range(len(df))

    plt.figure(figsize=(7, 0.35 * len(df) + 2), dpi=300)

    if se_metric is not None:
        # 95% CI = metric ± 1.96 * se
        xerr = 1.96 * df[se_metric]
    else:
        xerr = None

    plt.errorbar(df[metric], y, xerr=xerr, fmt='o', capsize=3)
    plt.yticks(y, df['study_id'])
    plt.axvline(df[metric].mean(), linestyle='--', label='Mean (unweighted)')
    plt.xlabel(xlabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


# -------------------------------------------------------------------
# 7. PRISMA flow diagram (simple)
# -------------------------------------------------------------------

def plot_prisma_flow(n_identified,
                     n_after_duplicates,
                     n_screened,
                     n_fulltext,
                     n_excluded_title_abs,
                     n_excluded_fulltext,
                     n_included,
                     filename: Path):
    """
    Dibuja un diagrama PRISMA muy simple con matplotlib.
    Rellena los números con tu flujo real.
    """
    plt.figure(figsize=(8, 10), dpi=300)
    ax = plt.gca()
    ax.axis('off')

    # Helper para dibujar un cuadro
    def box(x, y, w, h, text):
        rect = plt.Rectangle((x, y), w, h, fill=False)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', wrap=True)

    # Coordenadas (muy simples, puedes afinarlas)
    w, h = 3, 1.2

    # Cajas
    box(2.5, 8, w, h, f'Identified records\n(n={n_identified})')
    box(2.5, 6, w, h, f'Records after duplicates removed\n(n={n_after_duplicates})')
    box(2.5, 4, w, h, f'Titles/abstracts screened\n(n={n_screened})')
    box(0.2, 4, w, h, f'Records excluded\n(n={n_excluded_title_abs})')
    box(2.5, 2, w, h, f'Full-text articles assessed\n(n={n_fulltext})')
    box(0.2, 2, w, h, f'Full-text articles excluded\n(n={n_excluded_fulltext})')
    box(2.5, 0, w, h, f'Studies included in\nqualitative/quantitative synthesis\n(n={n_included})')

    # Flechas
    def arrow(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->"))

    arrow(2.5 + w/2, 8, 2.5 + w/2, 7.2)
    arrow(2.5 + w/2, 6, 2.5 + w/2, 5.2)
    arrow(2.5 + w/2, 4, 2.5 + w/2, 3.2)
    arrow(2.5 + w/2, 2, 2.5 + w/2, 1.2)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()


# -------------------------------------------------------------------
# 8. Exportar tablas estándar
# -------------------------------------------------------------------

def export_tables(df: pd.DataFrame,
                  summary_all: str,
                  summary_clinical: str,
                  summary_invitro: str):
    # Tabla 1: Características de estudios
    cols_char = ['study_id', 'setting', 'n_diseased', 'n_non_diseased',
                 'sensitivity', 'specificity', 'auc']
    cols_char = [c for c in cols_char if c in df.columns]
    df_char = df[cols_char].copy()
    df_char.to_csv(OUTPUT_DIR / "table1_study_characteristics.csv", index=False)

    # Tabla 2: Datos 2x2 + Se/Sp recalculadas
    cols_acc = ['study_id', 'TP', 'FN', 'FP', 'TN',
                'sens_recalc', 'spec_recalc']
    cols_acc = [c for c in cols_acc if c in df.columns]
    df_acc = df[cols_acc].copy()
    df_acc.to_csv(OUTPUT_DIR / "table2_accuracy_raw.csv", index=False)

    # Summary de modelos Reitsma en texto (puedes luego copiar a Word)
    with open(OUTPUT_DIR / "summary_reitsma_all.txt", "w") as f:
        f.write(summary_all)

    with open(OUTPUT_DIR / "summary_reitsma_clinical.txt", "w") as f:
        f.write(summary_clinical)

    with open(OUTPUT_DIR / "summary_reitsma_invitro.txt", "w") as f:
        f.write(summary_invitro)


# -------------------------------------------------------------------
# 9. Main
# -------------------------------------------------------------------

def main():
    # 1. Setup R + mada
    mada, base = setup_r_mada()
    print("Paquete 'mada' cargado en R.")

    # 2. Cargar datos y preparar 2x2
    df = load_and_prepare_data(DATA_PATH)
    print(f"{len(df)} estudios cargados.")

    # Separar clínicos / in vitro
    df_clinical = df[df['setting'] == 'clinical'].copy()
    df_invitro = df[df['setting'] == 'in_vitro'].copy()

    print(f"{len(df_clinical)} estudios clínicos, {len(df_invitro)} estudios in vitro.")

    # 3. Convertir a data.frame de R
    r_df_all = pandas2ri.py2rpy(df)
    r_df_clin = pandas2ri.py2rpy(df_clinical)
    r_df_inv = pandas2ri.py2rpy(df_invitro)

    # 4. Ejecutar Reitsma (global, sin moderadores)
    fit_all, summary_all = run_reitsma(mada, r_df_all, moderators=None)
    print("=== Reitsma: todos los estudios ===")
    print(summary_all)

    # 5. Reitsma por subgrupos (sin moderadores)
    fit_clin, summary_clinical = run_reitsma(mada, r_df_clin, moderators=None)
    fit_inv, summary_invitro = run_reitsma(mada, r_df_inv, moderators=None)

    # 6. Meta-regresión de ejemplo (setting +, por ejemplo, year)
    #    OJO: 'setting' es categórica; mada la tratará como factor.
    #    Asegúrate de que 'year' existe en el df si lo usas.
    moderators = []
    if 'setting' in df.columns:
        moderators.append('setting')
    if 'year' in df.columns:
        moderators.append('year')

    if moderators:
        print(f"Ejecutando meta-regresión con moderadores: {moderators}")
        fit_meta_reg, summary_meta_reg = run_reitsma(mada, r_df_all, moderators=moderators)
        with open(OUTPUT_DIR / "summary_reitsma_meta_regression.txt", "w") as f:
            f.write(summary_meta_reg)

    # 7. SROC global
    sroc_all = get_sroc_df(mada, fit_all)
    sroc_all.to_csv(OUTPUT_DIR / "sroc_points_all.csv", index=False)
    plot_sroc(sroc_all,
              title="SROC – All studies",
              filename=OUTPUT_DIR / "sroc_all.png")

    # 8. SROC clínico vs in vitro
    if len(df_clinical) > 0 and len(df_invitro) > 0:
        sroc_clin = get_sroc_df(mada, fit_clin)
        sroc_inv = get_sroc_df(mada, fit_inv)

        # Plot comparativo
        plt.figure(figsize=(6, 6), dpi=300)
        plt.plot(sroc_clin['fpr'], sroc_clin['sens'], label='Clinical SROC')
        plt.plot(sroc_inv['fpr'], sroc_inv['sens'], label='In vitro SROC')
        plt.plot([0, 1], [0, 1], linestyle='--', label='No discrimination')
        plt.xlabel('False Positive Rate (1 - Specificity)')
        plt.ylabel('Sensitivity')
        plt.title('SROC – Clinical vs In vitro')
        plt.legend()
        plt.grid(True)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / "sroc_clinical_vs_invitro.png", dpi=300)
        plt.close()

    # 9. Forest plots (Sens, Spec, AUC) – global y por subgrupos
    # Global
    forest_plot_metric(df, 'sens_recalc', 'se_sens',
                       title='Sensitivity – All studies',
                       xlabel='Sensitivity',
                       filename=OUTPUT_DIR / "forest_sensitivity_all.png")

    forest_plot_metric(df, 'spec_recalc', 'se_spec',
                       title='Specificity – All studies',
                       xlabel='Specificity',
                       filename=OUTPUT_DIR / "forest_specificity_all.png")

    if 'auc' in df.columns:
        forest_plot_metric(df, 'auc', None,
                           title='AUC – All studies',
                           xlabel='AUC',
                           filename=OUTPUT_DIR / "forest_auc_all.png")

    # Clínicos
    if len(df_clinical) > 0:
        forest_plot_metric(df_clinical, 'sens_recalc', 'se_sens',
                           title='Sensitivity – Clinical studies',
                           xlabel='Sensitivity',
                           filename=OUTPUT_DIR / "forest_sensitivity_clinical.png")

        forest_plot_metric(df_clinical, 'spec_recalc', 'se_spec',
                           title='Specificity – Clinical studies',
                           xlabel='Specificity',
                           filename=OUTPUT_DIR / "forest_specificity_clinical.png")

    # In vitro
    if len(df_invitro) > 0:
        forest_plot_metric(df_invitro, 'sens_recalc', 'se_sens',
                           title='Sensitivity – In vitro studies',
                           xlabel='Sensitivity',
                           filename=OUTPUT_DIR / "forest_sensitivity_invitro.png")

        forest_plot_metric(df_invitro, 'spec_recalc', 'se_spec',
                           title='Specificity – In vitro studies',
                           xlabel='Specificity',
                           filename=OUTPUT_DIR / "forest_specificity_invitro.png")

    # 10. PRISMA (rellena con tus números reales)
    plot_prisma_flow(
        n_identified=500,              # reemplaza con tus números
        n_after_duplicates=420,
        n_screened=420,
        n_fulltext=80,
        n_excluded_title_abs=340,
        n_excluded_fulltext=60,
        n_included=len(df),
        filename=OUTPUT_DIR / "prisma_flow.png"
    )

    # 11. Exportar tablas estándar
    export_tables(df, summary_all, summary_clinical, summary_invitro)

    print(f"Listo. Gráficos y tablas en: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
