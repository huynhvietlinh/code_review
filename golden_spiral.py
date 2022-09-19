def draw_golden_spiral (n):
  # every iteration turns 90 degrees
  angles = [0, 270, 180, 90]
  # this describes the direction 
  vectors = [(-1,0), (0,1), (1,0), (0,-1)]
  #
  f = [0 for i in range(n+1)]
  f[0] = 1
  f[1] = 1
  arc_center = (0,0)
  for i in range(1,n):
    plt.gca().add_patch(Arc(
      xy = arc_center,
      width = 2*f[i], height = 2*f[i],
      angle = angles[i%len(angles)],
      theta1 = 0, theta2 = 90, lw = 2    
    ))
    arc_center = (arc_center[0] + 
                 vectors[i%len(vectors)][0]*f[i-1], 
                 arc_center[1] + 
                 vectors[i%len(vectors)][1]*f[i-1])
    f[i+1] = f[i] + f[i-1]  
  plt.axis([-f[n], f[n], -f[n], f[n]])			
  plt.show()
