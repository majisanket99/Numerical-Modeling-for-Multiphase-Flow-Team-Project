import pygame


def screen_to_domain(mouse_x, mouse_y, window_width, window_height, domain_length_x, domain_length_y):
    """
    Convert screen pixel coordinates to simulation domain coordinates.

    Screen origin: top-left
    Simulation origin: bottom-left
    """
    x = mouse_x / window_width * domain_length_x
    y = domain_length_y - mouse_y / window_height * domain_length_y

    return x, y


def handle_events(ibm_object):
    """
    Handle keyboard and window events.

    Returns:
        running: False if the window should close
        reset_requested: True if user presses R
    """
    running = True
    reset_requested = False

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_ESCAPE:
                running = False

            if event.key == pygame.K_r:
                reset_requested = True

            if event.key == pygame.K_PLUS or event.key == pygame.K_EQUALS:
                ibm_object.increase_radius()

            if event.key == pygame.K_MINUS:
                ibm_object.decrease_radius()

    return running, reset_requested