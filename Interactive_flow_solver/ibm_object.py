import numpy as np


class IBMObject:
    """
    Circular immersed boundary object represented by Lagrangian boundary points.
    """

    def __init__(self, center_x, center_y, radius, n_points=80):
        self.center = np.array([center_x, center_y], dtype=float)
        self.previous_center = self.center.copy()

        self.radius = radius
        self.n_points = n_points

        self.boundary_points = self.compute_boundary_points()

    def compute_boundary_points(self):
        """
        Compute Lagrangian boundary points on a circular IBM object.
        """
        angles = np.linspace(0.0, 2.0 * np.pi, self.n_points, endpoint=False)

        x = self.center[0] + self.radius * np.cos(angles)
        y = self.center[1] + self.radius * np.sin(angles)

        return np.column_stack((x, y))

    def update_position(self, target_x, target_y, smoothing=0.15):
        """
        Move the object smoothly toward a target position.

        smoothing = 1.0 means instant movement.
        smoothing = 0.1 means slow smooth movement.
        """
        self.previous_center = self.center.copy()

        target = np.array([target_x, target_y], dtype=float)
        self.center = self.center + smoothing * (target - self.center)

        self.boundary_points = self.compute_boundary_points()

    def get_velocity(self, dt):
        """
        Compute object translational velocity from center movement.
        """
        if dt <= 0:
            return np.array([0.0, 0.0])

        return (self.center - self.previous_center) / dt

    def set_radius(self, new_radius):
        """
        Change radius and update boundary points.
        """
        self.radius = max(0.01, float(new_radius))
        self.boundary_points = self.compute_boundary_points()

    def increase_radius(self, amount=0.02):
        self.set_radius(self.radius + amount)

    def decrease_radius(self, amount=0.02):
        self.set_radius(self.radius - amount)