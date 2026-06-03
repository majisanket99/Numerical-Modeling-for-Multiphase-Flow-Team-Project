import pygame


def domain_to_screen(x, y, window_width, window_height, domain_length_x, domain_length_y):
    """
    Convert simulation domain coordinates to screen pixel coordinates.
    """
    sx = int(x / domain_length_x * window_width)
    sy = int((domain_length_y - y) / domain_length_y * window_height)

    return sx, sy


def draw_ibm_object(screen, ibm_object, window_width, window_height, domain_length_x, domain_length_y):
    """
    Draw IBM object center and Lagrangian boundary points.
    """
    center_x, center_y = ibm_object.center

    center_screen = domain_to_screen(
        center_x,
        center_y,
        window_width,
        window_height,
        domain_length_x,
        domain_length_y
    )

    # Draw transparent-looking filled circle approximation
    pygame.draw.circle(screen, (255, 180, 50), center_screen, 6)

    # Draw Lagrangian boundary points
    for point in ibm_object.boundary_points:
        px, py = point

        sx, sy = domain_to_screen(
            px,
            py,
            window_width,
            window_height,
            domain_length_x,
            domain_length_y
        )

        pygame.draw.circle(screen, (255, 220, 100), (sx, sy), 3)


def draw_text(screen, text, x, y, size=22, color=(230, 230, 230)):
    """
    Draw text on the Pygame screen.
    """
    font = pygame.font.SysFont(None, size)
    rendered_text = font.render(text, True, color)
    screen.blit(rendered_text, (x, y))