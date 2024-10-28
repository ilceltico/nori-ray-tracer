import mitsuba as mi
mi.set_variant("scalar_rgb")
scene = mi.load_file("/Users/stella/repos/cs440-2023-ilceltico/scenes/pa5/hacker3/cbox_mitsuba.xml")
image = mi.render(scene, spp=128)
mi.util.write_bitmap("/Users/stella/repos/cs440-2023-ilceltico/scenes/pa5/hacker3/cbox_mitsuba.exr", image)
