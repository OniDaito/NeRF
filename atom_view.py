import pyglet

config = pyglet.gl.Config(alpha_size=8)

window = pyglet.window.Window(config=config)

label = pyglet.text.Label('Hello, world',
  font_name='Times New Roman',
  font_size=36,
  x=window.width//2, y=window.height//2,
  anchor_x='center', anchor_y='center')

@window.event
def on_draw():
  window.clear()
  label.draw()

def launch():
  pyglet.app.run()

