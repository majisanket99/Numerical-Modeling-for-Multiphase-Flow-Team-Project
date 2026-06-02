import pygame
import numpy as np

# Window size
WIDTH, HEIGHT = 900, 500

# Simulation domain size
Lx, Ly = 4.0, 2.0

# Initialize pygame
pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Interactive IBM Flow Solver Test")

clock = pygame.time.Clock()

def screen_to_domain(mouse_x, mouse_y):
    """
    Convert mouse pixel coordinates to simulation domain coordinates.
    Screen origin is top-left.
    Simulation origin is bottom-left.
    """
    x = mouse_x / WIDTH * Lx
    y = Ly - mouse_y / HEIGHT * Ly
    return x, y

def domain_to_screen(x, y):
    """
    Convert simulation domain coordinates to screen pixel coordinates.
    """
    sx = int(x / Lx * WIDTH)
    sy = int((Ly - y) / Ly * HEIGHT)
    return sx, sy

running = True

while running:
    screen.fill((20, 20, 25))

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    # Read mouse position
    mouse_x, mouse_y = pygame.mouse.get_pos()
    object_x, object_y = screen_to_domain(mouse_x, mouse_y)

    # Draw IBM object following the mouse
    sx, sy = domain_to_screen(object_x, object_y)
    pygame.draw.circle(screen, (255, 180, 50), (sx, sy), 20)

    # Display text
    font = pygame.font.SysFont(None, 24)
    text = font.render(
        f"IBM object position: x={object_x:.2f}, y={object_y:.2f}",
        True,
        (230, 230, 230)
    )
    screen.blit(text, (20, 20))

    pygame.display.flip()
    clock.tick(60)

pygame.quit()