from draw_geom_funcs import *

def draw_geometry():
    # Target chamber
    draw_rectangle(-5.3, 5.3, -120.0, -55.0)
    draw_rectangle(-5.3, 5.3,  -55.0,  45.0)

    # High density concrete
    draw_rectangle(   5.3,  175.0, -120.0, -60.0)
    draw_hriz_line(-120.0, -125.0,   -5.3)
    draw_hriz_line( -60.0,  -65.0,   -5.3)
    draw_hriz_line(  60.0,  -65.0,  175.0)
    draw_hriz_line( 120.0, -125.0,  175.0)
    draw_vert_line(-125.0, -120.0,  120.0)
    draw_vert_line( -65.0,  -60.0,   60.0)
    draw_vert_line( 175.0, -120.0,  -60.0)
    draw_vert_line( 175.0,   60.0,  120.0)

    # Beam shaping assembly
    draw_rectangle(-40.0,  -5.3, -35.0, 35.0)
    draw_rectangle(  5.3,  40.0, -35.0, 35.0)
    draw_rectangle( 40.0,  40.3, -35.0, 35.0)
    draw_rectangle( 40.3,  45.3, -35.0, 35.0)
    draw_line     ( 45.3,  95.3, -35.0, -7.5)
    draw_line     ( 45.3,  95.3,  35.0,  7.5)
    draw_vert_line( 95.3,  -7.5, -40.0)
    draw_vert_line( 95.3,   7.5,  40.0)
    draw_hriz_line(-40.0,  95.3, 175.0)
    draw_hriz_line( 40.0,  95.3, 175.0)
    draw_vert_line(175.0, -40.0, -60.0)
    draw_vert_line(175.0,  40.0,  60.0)
