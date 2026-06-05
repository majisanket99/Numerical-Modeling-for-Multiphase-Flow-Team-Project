Project goal
Software requirements
How to build the Fortran solver
How to run the live interactive solver
How to run baseline and perturbed cases
How to generate comparison plots
Known limitations

This project adapts a Fortran IBM2D solver for real-time interactive flow control.
Mouse coordinates are read using xdotool and passed to the immersed boundary module.
The solver uses gnuplot X11 for live visualization.
The IBM object follows the cursor and introduces a localized disturbance into the flow.
Postprocessing scripts compare baseline and perturbed slice data.
