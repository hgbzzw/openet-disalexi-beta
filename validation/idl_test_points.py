def selector(x1, y1, x0=668565, y0=4636125, cs=30, rows=3000, cols=3000,
             tile_elements=255000):
    """Identify test point IDL tile and element from UTM coordinates

    The defaut parameters for this function hard hardcoded for the
    Mead validation dataset.
    """
    print('\n{} {}'.format(x1, y1))
    dx = abs(x1 - x0) / 30
    dy = abs(y1 - y0) / 30
    print('dx:       {}'.format(dx))
    print('dy:       {}'.format(dy))

    elements = dy * rows + dx
    print('element: {}'.format(elements))

    tile = float(elements) / tile_elements
    print('tile:    {}'.format(int(tile)))

    index = elements - int(tile) * tile_elements
    print('index:   {}'.format(int(index)))


if __name__ == '__main__':
    # AmeriFlux test points
    # Move test point to UL corner
    selector(x1=711690 - 15, y1=4560150 + 15)
    selector(x1=712260 - 15, y1=4560150 + 15)
    selector(x1=714750 - 15, y1=4561860 + 15)
