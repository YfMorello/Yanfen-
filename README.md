# Yanfen-
# Hi, I'm Yanfen
🎓 MSc Transport & Logistics | Data Analyst | Data & Analytics
💻 Python • Julia • R • CPLEX • SAS • SQL   
🌍 Passionate about using data to improve transport, logistics, and urban mobility, with an industrial background in shipping, logistics, supply chain management, and transport planning.  


### 🧠 Projects
- [Transport Data Explorer in Copenhagen](Urban-Transport-Optimisation) – Using operations research e.g. mixed-integer programming to optimise multimodal public transport and integrate them into a comprehensive solution for greener and more sustainable urban mobility.

## Career Goals
I aim to grow as a **Data Analyst** and transition into **Data Scientist** roles, applying analytical and data-driven decision-making to real-world industrial challenges in transport, logistics, and supply chain management.  

  
## 📫 Contact Me
Feel free to explore my projects and reach out with any questions or feedback. I'm excited to connect with like-minded professionals and recruiters in the data analytics field.  
- **LinkedIn:** https://www.linkedin.com/in/yanfen-chen-a0a772252
- **Email:** yanfen.chen08@gmail.com


# Integrating Bike-Sharing Systems (BSS) with Bus Rapid Transit (BRT) Along the Ring 4 Corridor (O4)

A data-driven operations research and network optimization project focused on multimodal urban mobility planning in the Greater Copenhagen area. This repository contains the complete codebase, data preprocessing pipelines, and optimization models developed to evaluate and optimize the integration of bike-sharing dock locations with the Ring 4 Bus Rapid Transit (BRT) corridor (specifically targeting the 400S bus line and surrounding mobility hubs).

---

## 📌 Project Overview

Urban public transit often suffers from the "first-and-last-mile" problem. This project investigates how strategic placement of bike-sharing stations (BSS) along high-capacity transit corridors like Copenhagen's Ring 4 can enhance feeder transit efficiency, reduce total travel time, and optimize capital allocation for urban planners.

### Key Objectives:
* **Multimodal Data Integration**: Processed large-scale Origin-Destination (OD) trip matrices, spatial GTFS/bus stop coordinates, and road network distances.
* **Mathematical Modeling & Optimization**: Implemented optimization formulations (including facility location and covering models) in Julia to determine optimal bike dock placement under varying budget/capacity parameters ($M$).
* **Sensitivity & Parameter Analysis**: Conducted comparative analyses across different scenario baselines to evaluate variable interactions between bike-sharing infrastructure and bus corridor demand.

---

## 🏗️ Repository Architecture

```text
├── data/                    # Open-source mobility data, bus stops, and processed OD matrices
├── src/
│   ├── preprocessing.py     # Python scripts for data cleaning, spatial filtering, and OD shortest paths
│   └── optimization.jl      # Julia code for mathematical formulation and location-allocation optimization
├── outputs/                 # Generated maps, distribution plots, and sensitivity analysis results
├── requirements.txt         # Python dependencies
└── README.md                # Project documentation
