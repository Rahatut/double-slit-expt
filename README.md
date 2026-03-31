# Young’s Double Slit Experiment Simulation

This project simulates the **Young’s double-slit experiment** in C++ using **Raylib**, visualizing electron interference patterns based on slit and screen distances. The simulation models the probabilistic behavior of electrons and demonstrates how interference patterns change with varying experimental parameters.

---

## Features

- Generates a **wavefunction** based on slit and screen distances.
- Calculates **intensity probabilities** over a range of firing angles.
- Normalizes probabilities using the **trapezoidal rule**.
- Fires a beam of electrons with **randomized positions** according to probability distributions.
- Plots electron impacts on a screen in real-time.
- Dynamically updates **slit distance and screen distance** to observe changes in the interference pattern.

---

## How the Simulation Works

1. **Initialize the experiment:** Create a new wavefunction with given slit and screen distances.  
2. **Wavefunction generation:** Calculate intensities over linearly spaced firing angles, storing them as probabilities.  
3. **Normalization:** Probabilities are normalized using the trapezoidal rule.  
4. **Electron firing:** A beam of `n` electrons is simulated. Each electron’s position `(x, y)` is randomly chosen based on the probability distribution.  
5. **Screen plotting:** Each electron is plotted as a circle on the screen, centered at the screen’s midpoint.  
6. **Dynamic update:** Slit distance and screen distance are gradually increased to show their effect on the interference pattern.

---

## Technologies Used

- **C++** – core programming language  
- **Raylib** – graphics library for visualization  
- **Git** – version control  

---

## Demo

A short video demonstrating the effect of increasing slit distance on the interference pattern is available [here](https://www.linkedin.com/posts/rahatut-tahrim_replicating-youngs-double-slit-experiment-activity-7253357501750878208-_eAb?utm_source=share&utm_medium=member_desktop&rcm=ACoAAEYnkQkBm-f0BsRx6jgpOBP-8e2drArgSn8).

---

## How to Run

1. Clone the repository:  
   ```bash
   git clone https://github.com/Rahatut/double-slit-expt.git
