def slope(h, w):
    #return (8*pi*G*rho*(w-1)-h[0]**2-h[0]*h[1]-h[0]*h[2])
    return ((w-1)-h[0]**2-h[0]*h[1]-h[0]*h[2])

def euler(step_size, step_count, initial_values, w, start = 0):
    vals = {start:initial_values}
    last = initial_values
    for n in range(step_count):
        hcs = [last[0]+slope(last, w)*step_size, last[1]+slope(last[1:]+[last[0]], w)*step_size,
                last[2]+slope([last[2]]+last[:2], w)*step_size]
        vals[round(start+(n+1)*step_size, 3)] = hcs
        last = hcs
    return vals

if __name__ == '__main__':
    import matplotlib.pyplot as p
    import pylab
    h = euler(.001, 100, [20, 20, 10], 0, .001)
    x = [i/float(1000) for i in range(1, 100)]
    y1 = [h[t][0] for t in x]
    y2 = [h[t][1] for t in x]
    y3 = [h[t][2] for t in x]
    p.scatter(x, y1, c = 'b')
    p.scatter(x, y2, c = 'g')
    p.scatter(x, y3, c = 'r')
    p.show()
