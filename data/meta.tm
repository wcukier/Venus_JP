KPL/MK

   This is the meta-kernel used for the Spring 2021 JP on
   Microbial Transfer from Venus by Wolf Cukier

   The names and contents of the kernels referenced by this
   meta-kernel are as follows:

   File name                   Contents
   --------------------------  -----------------------------
   naif0012.tls                Generic LSK
   de430.bsp                   Lunar/Planetary Ephemeris SCK



   \begindata
   PATH_VALUES = ( 'data/SPICE' )
   
   PATH_SYMBOLS = ( 'PTH')
   
   KERNELS_TO_LOAD = ( '$PTH/naif0012.tls',
                       '$PTH/de440s.bsp',
)
   \begintext