================================================================================
COMMENTS REGARDING THE DIFFERENT FITTING FUNCTIONS:

Find line => 'g'/'r'
Metallic line => ('g' for speed) 'r' -> 'vr_Z' -> 'vrg_Z'
Hydrogen/Helium line => ('l' for speed) 'v' -> 'vr_H' -> 'vrg_H'

Note: functions with lorentzians are defined with a variable "y" value for a
      better fitting with the continuum.

Pure gaussian profile ----------------------------------------------------------
  - Name  -> f_gaussian1(x,A,x0,sigma)
  - Alias -> 'g'
  - Lines -> + Z lines | - H/He lines
  - Range vsini/FWHM -> ~<60????/????A

Pure lorentzian profile --------------------------------------------------------
  - Name  -> f_lorentzian(x,A,x0,gamma)
  - Alias -> 'l'
  - Lines -> - Z lines | + H/He lines

  - Note: Very good for H/He lines with moderate rotation 100-250kms

Voigt profile (g x l) ----------------------------------------------------------
  - Name  -> f_voigt(x,A,x0,sigma,gamma,y)
  - Alias -> 'v'
  - Lines -> + Z lines | - H/He lines
  - Range vsini/FWHM -> ~< 150kms/3.8A for R ~   2500
                        ~<  60kms/2.0A for R ~< 15000
                       (~<  30kms/1.0A for R  > 15000)
  - Note: Ok for first rough estimation of the FWHM

Rotational profile (g x r) -----------------------------------------------------
  - Name  -> f_rot(x,A,x0,sigma,vsini)
  - Alias -> 'r'
  - Lines -> +++ Z lines | + H/He lines
  - Range vsini/FWHM -> ~<410kms/10A for any R
  - Note: Very good for metallic lines at any rotation
          Useful for first rough estimation of the FWHM

Voigt with rotation profile (g x l x r) ----------------------------------------
  - Name  -> f_voigtrot(x,A,x0,sigma,gamma,vsini,y)
UNIFICAR

  - Alias -> 'vr_Z' // 'vr_H'
  - Lines -> ++ Z lines | - H/He lines // - Z lines | ++ H/He lines
  - Range vsini/FWHM -> ~<160kms/4A for any R // ~>160kms/4A for any R
  - Note: Good for metallic lines with low-moderate rotation //
          Good for H/He lines with moderate-high rotation

Voigt with rotation profile plus gaussian (g x l x r + g) ----------------------
  - Name  -> f_vrg(x,A,x0,sigma,gamma,vsini,y)
  - Alias -> 'vrg_Z' // 'vrg_H'
  - Lines -> ++ Z lines | - H/He lines // - Z lines | +++ H/He lines
  - Range vsini/FWHM -> ~<410kms/10A for any R // ~>410kms/15A for any R
  - Note: Good for metallic lines with any rotation //
          Very good for H/He lines with any rotation
  - Note2: 'vrg_Z' should only be used if the line is complex and somehow
           improves 'vr/r'. The lower limit for 'y' parameter is set to -0.1
           and therefore could produce larger EWs then other functions.


================================================================================
FITLINE DIAGRAM:

+-----------------+
|Basic input data |
+--+--------------+
   |
+--v--------------+
|    Iteration    <----------------------^-----------------+
+-----------+-----+                      |                 |
            |                            |                 |
     +------v------------------+         |                 |
     | Last iteration reached? |         |                 |
     +---+----+----------------+         |                 |
         |    |                          |                 |
         v    v                          |                 |
+------+YES  NOP                         |                 |
|             +                          |                 |
|             |                          |                 |
|   +---------v---+                      |                 |
|   | Resampling? |                      |                 |
|   +--+------+---+                      |                 |
|      |      |                          |                 |
|      v      v    +----------------+    |                 |
|     NOP    YES+--> Resampled spec |    |                 |
|      +           +-------+--------+    |                 |
|      |                   |             |                 |
|      +<------------------+             |                 |
|      |                                 |                 |
|   +--v-------------------+             |                 |
|   | Find normalization   |             |                 |
|   | regions x3 times     |             |                 |
|   +--+-------------------+             |                 |
|      |                                 |                 |
|   +--v-------------------+             |                 |
|   | Normalizing the flux |             |                 |
|   +--+-------------------+             |                 |
|      |                          +------+--------------+  |
|   +--v-------------------+      | Sets new best width |  |
|   | Try the line fitting |      | and fitting results |  |
|   +--+----------+--------+      +-----------^---------+  |
|      |          |                           |            |
|      v          v         +------------+    |            |
|     BAD       GOOD+-------> Calculates |    +            |
|      +                    | the FWHM   |   YES           |
|      |                    +-----+------+    ^            |
|   +--v---------------+          |           |            |
|   | First iteration? <----+  +--v-----------+----------+ |
|   +--+----------+----+    |  |                         | |
|      |          |         |  | Resolution < FWHM < 15? | |
|      v          v         |  |                         | |
|     YES        NOP        |  +--+----------------------+ |
|      |          +         |     |                        |
|  +---v---+      |         |     v                        |
|  | BREAK |      |         +---+NOP                       |
|  +-------+      |                                        |
|                 |                                        |
|  +--------------+-------+                                |
|  | Reach last iteration +--------------------------------+
|  +----------------------+
|
|  +---------------------+
+--> Finding line center |
   +--+------------------+
      |
   +--v---------------------+
   | Line within tolerance? |
   +---+----------+---------+
       |          |
       v          v
      NOP        YES
       |          +      +--------------------------------+
   +---v---+      |      | Calculates the EW, FWHM, SNR   |
   | BREAK |      +----->+ and the quality of the fitting |
   +-------+             +--------------------------------+

================================================================================
