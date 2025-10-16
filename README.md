# 🔥 ThermalFVT_Homogenization

Este repositório fornece implementações em **MATLAB** voltadas à **homogeneização térmica de materiais heterogêneos**, com base na **Teoria dos Volumes Finitos (FVT)**.  
O objetivo é determinar a **matriz de condutividade térmica efetiva** de um meio periódico, permitindo comparar diferentes formulações numéricas e analisar a resposta térmica global a partir do comportamento microscópico.

Duas formulações distintas são apresentadas:

- **Energy-based formulation** — baseada na equivalência energética;
- **Mean-field formulation** — beseada teoria de campos médios.

### 🧩 Recursos principais

- Periodic boundary conditions automatically enforced.
- Homogenization of periodic cellular materials based on the concept of Repeating Unit Cell (RUC).
- Computation of fluctuating and surface-averaged temperatures on the subvolume faces.
- Computation of the effective thermal conductivity matrix (**k***).

## ⚙️ Getting Started

Save the program **ThermalFVT_energy_based.m** (or **ThermalFVT_mean_field.m**) and start MATLAB in the same directory.  
The program can be executed with the following command:

```matlab
ThermalFVT_energy_based(nx, ny, frac, k_m, k_i)
```
where:
- **nx** and **ny** define the number of subvolumes along the x and y directions, respectively;
- **frac** is the volume fraction of the inclusion;
- **k_m** is the thermal conductivity of the matrix;
- **k_i** is the thermal conductivity of the inclusion.

Run the main function as shown below:
```matlab
ThermalFVT_energy_based(200, 200, 0.2, 1, 50)   % Energy-based formulation
ThermalFVT_mean_field(200, 200, 0.2, 1, 50)     % Mean-field formulation
```
## 👨‍💻 Authors

- **Diogo Tiago dos Santos** — [diogo.santos@ctec.ufal.br](mailto:diogo.santos@ctec.ufal.br)  
- **Márcio André Araújo Cavalcante** — [marcio.cavalcante@ceca.ufal.br](mailto:marcio.cavalcante@ceca.ufal.br)  
- **Romildo dos Santos Escarpini Filho** — [romildo.escarpini@penedo.ufal.br](mailto:romildo.escarpini@penedo.ufal.br)  
- **Arnaldo dos Santos Júnior** — [arnaldo@ctec.ufal.br](mailto:arnaldo@ctec.ufal.br)
