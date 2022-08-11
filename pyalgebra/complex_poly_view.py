import pyglet
import moderngl
import numpy as np
import random
import math



class Renderer():    
    def render(self, width, height, center, scale):
        raise NotImplementedError



class ComplexFunc(Renderer):
    def __init__(self, ctx, coeffs):
        super().__init__()

        assert len(coeffs) >= 2
        
        self.prog = ctx.program(
            vertex_shader = """
                #version 430

                in vec2 in_vert;
                out vec2 v_coord;

                uniform vec2 size;
                uniform vec2 center;
                uniform float scale;

                void main() {
                    float side = sqrt(size.x * size.y);
                    gl_Position = vec4(in_vert, 0.0, 1.0);
                    v_coord = scale * size * in_vert / side - center;
                }
            """,
            fragment_shader = """
                #version 430

                in vec2 v_coord;
                out vec4 f_colour;
                uniform float coeffs_r[""" + str(len(coeffs)) + """];
                uniform float coeffs_i[""" + str(len(coeffs)) + """];
                int degree = """ + str(len(coeffs) - 1) + """;

                vec3 hsl2rgb( in vec3 c )
                {
                    vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );

                    return c.z + c.y * (rgb-0.5)*(1.0-abs(2.0*c.z-1.0));
                }

                vec2 product(vec2 z, vec2 w) {
                    return vec2(z.x * w.x - z.y * w.y, z.x * w.y + z.y * w.x);
                }

                void main() {
                    vec2 z = vec2(coeffs_r[degree], coeffs_i[degree]);

                    for (int i = degree - 1; i >= 0; i --){
                        z = vec2(coeffs_r[i], coeffs_i[i]) + product(v_coord, z);
                    }

                    float tau = 6.283185307179586;
                    f_colour = vec4(hsl2rgb(vec3(atan(z.y, z.x) / tau, 1, 0.5 + 0.1 * sin(3 * log(length(z))))), 1);
                }
            """
        )
        

        coeffs = [complex(c) for c in coeffs]

        self.prog["coeffs_r"].value = [c.real for c in coeffs]
        self.prog["coeffs_i"].value = [c.imag for c in coeffs]

        vertices = np.array([-1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0])

        self.vbo = ctx.buffer(vertices.astype('f4'))
        self.vao = ctx.simple_vertex_array(self.prog, self.vbo, 'in_vert')

        
    def render(self, width, height, center, scale):
        self.prog['center'].value = center
        self.prog['scale'].value = scale
        self.prog['size'].value = (width, height)

        self.vao.render(moderngl.TRIANGLE_STRIP)



        

class ImplicitFunc(Renderer):
    def __init__(self, ctx, glsl_func):
        super().__init__()
        
        self.prog = ctx.program(
            vertex_shader = """
                #version 430

                in vec2 in_vert;
                out vec2 v_coord;

                uniform vec2 size;
                uniform vec2 center;
                uniform float scale;

                void main() {
                    float side = sqrt(size.x * size.y);
                    gl_Position = vec4(in_vert, 0.0, 1.0);
                    v_coord = scale * size * in_vert / side - center;
                }
            """,
            fragment_shader = """
                #version 430

                in vec2 v_coord;
                uniform float radius;
                out vec4 f_colour;

                vec3 hsl2rgb( in vec3 c )
                {
                    vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );
                    return c.z + c.y * (rgb-0.5)*(1.0-abs(2.0*c.z-1.0));
                }

                vec2 product(vec2 z, vec2 w) {
                    return vec2(z.x * w.x - z.y * w.y, z.x * w.y + z.y * w.x);
                }

                float func(float x, float y) {
                    return """ + glsl_func + """;
                }

                void main() {
                    float eps1 = radius;
                    float eps2 = radius / 1.41421356237;
                    //float y = func(v_coord);

                    bool a = func(v_coord.x + eps1, v_coord.y) < 0;
                    bool b = func(v_coord.x - eps1, v_coord.y) < 0;
                    bool c = func(v_coord.x, v_coord.y + eps1) < 0;
                    bool d = func(v_coord.x, v_coord.y - eps1) < 0;
                    bool e = func(v_coord.x + eps2, v_coord.y + eps2) < 0;
                    bool f = func(v_coord.x + eps2, v_coord.y - eps2) < 0;
                    bool g = func(v_coord.x - eps2, v_coord.y + eps2) < 0;
                    bool h = func(v_coord.x - eps2, v_coord.y - eps2) < 0;

                    if ((a && b && c && d && e && f && g && h) || (!a && !b && !c && !d && !e && !f && !g && !h)) {
                        discard;
                    } else {
                        f_colour = vec4(0, 0, 0, 1);
                    }
                }
            """
        )

        vertices = np.array([-1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0])

        self.vbo = ctx.buffer(vertices.astype('f4'))
        self.vao = ctx.simple_vertex_array(self.prog, self.vbo, 'in_vert')

        
    def render(self, width, height, center, scale):
        pixel_rad = 3
        
        self.prog['center'].value = center
        self.prog['scale'].value = scale
        self.prog['radius'].value = pixel_rad * scale / math.sqrt(width * height)
        self.prog['size'].value = (width, height)

        self.vao.render(moderngl.TRIANGLE_STRIP)



class ZoomWindow(pyglet.window.Window):
    def __init__(self, size, renderers):
        super().__init__(*size, resizable = True)
        self.ctx = moderngl.create_context()
        self.ctx.enable_only(moderngl.BLEND)

        renderers = [renderer(self.ctx) for renderer in renderers]
        
        for renderer in renderers:
            assert issubclass(type(renderer), Renderer)

        self.scale = 1.0
        self.center = (0.0, 0.0)
        self.renderers = renderers

    def convert_inwards(self, x, y):
        side = math.sqrt(self.width * self.height)
        return self.center[0] + self.width * self.scale * (2 * x / self.width - 1) / side, self.center[1] + self.height * self.scale * (2 * y / self.height - 1) / side
    def convert_outwards(self, x, y):
        raise NotImplementedError

    def on_draw(self):
        self.clear()
        self.render()

    def on_mouse_scroll(self, x, y, dx, dy):
        p1 = self.convert_inwards(x, y)
        self.scale *= 0.9 ** dy
        p2 = self.convert_inwards(x, y)
        self.center = tuple(self.center[i] - p1[i] + p2[i] for i in [0, 1])

    def render(self):
        self.ctx.clear()
        for renderer in self.renderers:
            renderer.render(self.width, self.height, self.center, self.scale)








class Fractal(Renderer):
    def __init__(self, ctx):
        super().__init__()
        self.ctx = ctx
        self.prog = self.ctx.program(
            vertex_shader = """
                #version 430

                in vec2 in_vert;
                out vec2 v_text;

                uniform vec2 size;

                void main() {
                    float side = sqrt(size.x * size.y);
                    gl_Position = vec4(in_vert, 0.0, 1.0);
                    v_text = size * in_vert / side;
                }
            """,
            fragment_shader = """
                #version 430

                in vec2 v_text;
                out vec4 f_color;

                uniform vec2 center;
                uniform float scale;
                uniform int iter;

                vec3 hsl2rgb( in vec3 c )
                {
                    vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );

                    return c.z + c.y * (rgb-0.5)*(1.0-abs(2.0*c.z-1.0));
                }


                void main() {
                    vec2 c;
                    int i;

                    c.x = v_text.x * scale - center.x;
                    c.y = v_text.y * scale - center.y;

                    vec2 z = c;

                    for (i = 0; i < iter; i++) {
                        float x = (z.x * z.x - z.y * z.y) + c.x;
                        float y = (z.y * z.x + z.x * z.y) + c.y;

                        if ((x * x + y * y) > 4.0) {
                            break;
                        }

                        z.x = x;
                        z.y = y;
                    }

                    if (i == iter) {
                        discard;
                        f_color = vec4(0, 0, 0, 1);
                    } else {
                        f_color = vec4(hsl2rgb(vec3(i / 10.0, 1, 0.5)), 1);
                    }
                }
            """
        )

        vertices = np.array([-1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0])

        self.vbo = self.ctx.buffer(vertices.astype('f4'))
        self.vao = self.ctx.simple_vertex_array(self.prog, self.vbo, 'in_vert')

        
    def render(self, width, height, center, scale):
        self.prog['center'].value = center
        self.prog['iter'].value = 12
        self.prog['scale'].value = scale
        self.prog['size'].value = (width, height)

        self.vao.render(moderngl.TRIANGLE_STRIP)
        




def run(coeffs, rects):
    
    def implicit_rect(ax, bx, ay, by):
        ax, bx, ay, by = float(ax), float(bx), float(ay), float(by)
        mx = 0.5 * (bx + ax)
        my = 0.5 * (by + ay)
        EPS = 10 ** -10
        if abs(ax - bx) < EPS:
            ax -= EPS / 2
            bx += EPS / 2
        if abs(ay - by) < EPS:
            ay -= EPS / 2
            by += EPS / 2
        return lambda ctx : ImplicitFunc(ctx, f"max(abs(x - ({mx})) * {2 / (bx - ax)}, abs(y - ({my})) * {2 / (by - ay)}) - 1")

    window = ZoomWindow((1600, 1000), [lambda ctx : ComplexFunc(ctx, coeffs)] + [implicit_rect(*rect) for rect in rects])


    pyglet.app.run()


if __name__ == "__main__":
    run([-1, 0, 0, 1], [[-2, 2, -1, 1]])





















