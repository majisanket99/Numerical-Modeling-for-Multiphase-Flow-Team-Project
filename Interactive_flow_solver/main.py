import pygame

from ibm_object import IBMObject
from interaction import screen_to_domain, handle_events
from visualization import draw_ibm_object, draw_text


# Window size in pixels
WIDTH, HEIGHT = 900, 500

# Simulation domain size
Lx, Ly = 4.0, 2.0

# Time step for object velocity calculation
dt = 1.0 / 60.0


def main():
    pygame.init()

    screen = pygame.display.set_mode((WIDTH, HEIGHT))
    pygame.display.set_caption("Interactive IBM Flow Solver")

    clock = pygame.time.Clock()

    # Create IBM object at domain center
    ibm_object = IBMObject(
        center_x=Lx / 2.0,
        center_y=Ly / 2.0,
        radius=0.15,
        n_points=80
    )

    running = True

    while running:
        screen.fill((20, 20, 25))

        running, reset_requested = handle_events(ibm_object)

        if reset_requested:
            ibm_object = IBMObject(
                center_x=Lx / 2.0,
                center_y=Ly / 2.0,
                radius=0.15,
                n_points=80
            )

        # Mouse position in screen coordinates
        mouse_x, mouse_y = pygame.mouse.get_pos()

        # Convert mouse position to simulation domain
        target_x, target_y = screen_to_domain(
            mouse_x,
            mouse_y,
            WIDTH,
            HEIGHT,
            Lx,
            Ly
        )

        # Smoothly move IBM object toward mouse
        ibm_object.update_position(target_x, target_y, smoothing=0.15)

        # Compute object velocity
        object_velocity = ibm_object.get_velocity(dt)

        # Draw IBM object
        draw_ibm_object(screen, ibm_object, WIDTH, HEIGHT, Lx, Ly)

        # Draw information
        draw_text(
            screen,
            f"IBM center: x={ibm_object.center[0]:.2f}, y={ibm_object.center[1]:.2f}",
            20,
            20
        )

        draw_text(
            screen,
            f"Velocity: ux={object_velocity[0]:.2f}, uy={object_velocity[1]:.2f}",
            20,
            45
        )

        draw_text(
            screen,
            f"Radius: {ibm_object.radius:.2f} | Boundary points: {ibm_object.n_points}",
            20,
            70
        )

        draw_text(
            screen,
            "Controls: Move mouse | +/- change radius | R reset | ESC quit",
            20,
            HEIGHT - 35,
            size=20,
            color=(180, 180, 180)
        )

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()


if __name__ == "__main__":
    main()