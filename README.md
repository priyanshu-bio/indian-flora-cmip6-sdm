# 🧬 Predicting Climate-Driven Range Shifts of Key Indian Botanical Guilds

An ecological niche modeling framework using **Maxnet** and **CMIP6 climate scenarios** (SSP2-4.5 & SSP5-8.5) to evaluate habitat suitabilities for 10 ecologically and economically important Indian plant species.

## 🌿 Guilds & Target Species
- **Invasive Alien Species:** *Lantana camara*, *Parthenium hysterophorus*
- **Himalayan High-Altitude Medicinals:** *Aconitum heterophyllum*, *Picrorhiza kurroa*, *Rheum emodi*
- **Endemic/Native Medicinal Flora:** *Withania somnifera*, *Pterocarpus santalinus*, *Garcinia indica*, *Gymnema sylvestre*, *Tinospora cordifolia*

## 📊 Model Performance
- **Mean AUC:** 0.943
- **Mean TSS:** 0.786

## Sample Maps

<img width="2200" height="2200" alt="Gymnema_sylvestre_current_suitability_clean" src="https://github.com/user-attachments/assets/edac22d2-c75d-47b7-ad5e-d21d9eb6c009" />
<img width="2200" height="2200" alt="Parthenium_hysterophorus_current_suitability_clean" src="https://github.com/user-attachments/assets/505cdb54-bad4-41a1-b33e-758b02ae2763" />
<img width="2200" height="2200" alt="Pterocarpus_santalinus_current_suitability_clean" src="https://github.com/user-attachments/assets/a25d6adf-a4c4-4f4b-b7ed-86a4dacf9893" />


## 📁 Repository Structure
- `scripts/`: R scripts for preprocessing, Maxnet modeling, and future projection analysis.
- `outputs/`: Summary performance metrics and projected range-shift visualizations.

## 🚀 How to Run
1. Clone this repository:
   ```bash
   git clone [https://github.com/your-username/species-distribution-modeling-cmip6.git](https://github.com/your-username/species-distribution-modeling-cmip6.git)
   ```
2. Open the project in RStudio and install dependencies (`maxnet`, `terra`, `predicts`, `ggplot2`).
3. Execute scripts in sequence from `scripts/01_data_preprocessing.R`.

## 📜 Citation & License
This repository is released under the MIT License. 
If you use this code or workflow, please cite our associated manuscript under review at *Biodiversity and Conservation*.
