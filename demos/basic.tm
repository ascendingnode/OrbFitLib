   \begintext
   
     This is a basic metakernel for the solar system.

     Edit the string in PATH_VALUES to match where you store your kernels.

   \begindata

     PATH_VALUES  = ('../../../kernels/')

     PATH_SYMBOLS = ('KDIR')

     KERNELS_TO_LOAD = ( '$KDIR/naif0010.tls',
                         '$KDIR/de432s.bsp' )

   \begintext
