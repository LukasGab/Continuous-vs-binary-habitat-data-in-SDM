# Continuous-vs-binary-habitat-data-in-SDM
This repository was created as a supporting material for the article entitled: *"Habitats as predictors in species distribution models: Shall we use continuous or binary data?"*

Code authors: [Lukas Gabor](https://scholar.google.cz/citations?user=pLQXY5wAAAAJ&hl=cs) & [Petr Keil](https://scholar.google.cz/citations?user=SUAqa68AAAAJ&hl=cs&oi=ao)

Date the repository was created: 16/02/2022

**Relevant paper:**
[Gábor, L., Šímová, P., Keil, P., Zarzo‐Arias, A., Marsh, C. J., Rocchini, D., ... & Moudrý, V. (2022). Habitats as predictors in species distribution models: Shall we use continuous or binary data?. Ecography, e06022.](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.06022)

**Files Description:**

birds_name.csv: File with species Latin and English names. Used for results save and plot results. 

data.csv: File with species data and environmetal variables used in our study. Species data were obtain from the [Atlas of breeding bird distribution in the Czech Republic 2001-2003](https://scholar.google.cz/scholar?hl=cs&as_sdt=0%2C7&q=Atlas+of+breeding+bird+distribution+in+the+Czech+Republic+2001-2003&btnG=). Environmetal variables were downloaded from CORINE Land Cover [webpages](https://land.copernicus.eu/pan-european/corine-land-cover).

R_Script.R: Step-by-step guideline for loading data, create binary habitat variables adn run GLM models.
