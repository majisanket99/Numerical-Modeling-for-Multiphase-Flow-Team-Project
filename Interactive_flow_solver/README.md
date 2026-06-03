# Numerical-Modeling-for-Multiphase-Flow-Team-Project
Topic 4: Interactive real-time flow solver

# IBM Real-Time Flow Solver

This project develops an interactive real-time 2D flow solver using Python, NumPy, and Pygame.  
The goal is to control an immersed boundary object with the mouse and use it to perturb a vortex street.

## Setup

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python main.py

1. Mouse-controlled circle
2. Convert circle into IBM boundary points
3. Create a simple 2D grid using NumPy
4. Visualize a scalar field, for example vorticity
5. Add a simple flow field
6. Add IBM forcing around the circle
7. Move the IBM object using the mouse
8. Save time, object position, and maximum vorticity
9. Compare unperturbed and perturbed vortex street cases