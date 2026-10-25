# Data Directory / Directorio de Datos

**Important Note:** Due to GitHub's file size limits (100MB), the raw `.dta` files from the ENEMDU survey are **not** stored in this repository. 
*(Nota importante: Debido a los límites de tamaño de GitHub, los archivos originales .dta de la ENEMDU no se almacenan en este repositorio).*

## How to reproduce the analysis? / ¿Cómo reproducir el análisis?

1. Download the raw data from the official INEC portal: [https://www.ecuadorencifras.gob.ec/empleo-anual/](https://www.ecuadorencifras.gob.ec/empleo-anual/)
2. Create a folder named `raw` inside this `data` directory and place the downloaded files there (`temp_enemdu_2021.dta`, etc.).
3. The R script will automatically generate a `processed` folder with the clean data.
4. **The aggregated dataset used for the functional clustering is already available** in `outputs/tables/aggregated_data_46curves.csv`.
