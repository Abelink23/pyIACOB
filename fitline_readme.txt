================================================================================
COMMENTS REGARDING THE DIFFERENT FITTING FUNCTIONS:

Usage:
 To find line => 'g'/'r'
 To fit metallic line => ('g' for speed) 'r' -> 'vr_Z' -> 'vrg_Z'*
 To fit Hydrogen/Helium line => ('l' for speed) 'v' -> 'vr_H' -> 'vrg_H'

Note: functions with lorentzians are defined with a variable "y" value for a
      better fitting with the continuum.

Pure gaussian profile ----------------------------------------------------------
  - Name  -> f_gaussian1(x,A,x0,sigma)
  - Alias -> 'g'
  - Lines -> + Z lines | - H/He lines
  - Range vsini/FWHM -> ~<60kms/2.5A for R ~< 20000
  - Note: Good for first rough estimation of the FWHM in metallic lines

Pure lorentzian profile --------------------------------------------------------
  - Name  -> f_lorentzian(x,A,x0,gamma)
  - Alias -> 'l'
  - Lines -> - Z lines | + H/He lines
  - Range vsini/FWHM -> (~>100kms/3.5 for any R) - depends on the range
  - Note: Ok for first rough estimation of the FWHM in metallic lines
          Good for H/He lines with moderate rotation 100-200kms
          Good for H/He at lgf<2 + <100kms + R<5000 or lgf>2 + >250kms

Voigt profile (g x l) ----------------------------------------------------------
  - Name  -> f_voigt(x,A,x0,sigma,gamma,y)
  - Alias -> 'v'
  - Lines -> - Z lines | ++ H/He lines
  - Range vsini/FWHM -> ~>100kms/3.5 for any R
  - Note: Good for first rough estimation of the FWHM in H/He lines
          Very good for lgf<2 + >100kms
          For metallic lines better use the 'r' profile.

Rotational profile (g x r) -----------------------------------------------------
  - Name  -> f_rot(x,A,x0,sigma,vsini)
  - Alias -> 'r'
  - Lines -> +++ Z lines | - H/He lines
  - Range vsini/FWHM -> ~<410kms/10A for any R
  - Note: Very good for first rough estimation of the FWHM in metallic lines
          Very good for metallic lines at any rotation

Voigt with rotation profile (g x l x r) ----------------------------------------
  - Name  -> f_voigtrot(x,A,x0,sigma,gamma,vsini,y)
  - Alias -> 'vr_Z' // 'vr_H'
  - Lines -> ++ Z lines | - H/He lines // - Z lines | ++ H/He lines
  - Range vsini/FWHM -> ~<160kms/4A for any R // ~>160kms/4.5A for any R
  - Note: Good for metallic lines with low-moderate rotation //
          Good for H/He lines with moderate-high rotation

Voigt with rotation profile plus gaussian (g x l x r + g) ----------------------
  - Name  -> f_vrg(x,A,x0,sigma,gamma,vsini,y)
  - Alias -> 'vrg_Z' // 'vrg_H'
  - Lines -> ++ Z lines | - H/He lines // - Z lines | +++ H/He lines
  - Range vsini/FWHM -> ~<410kms/10A for any R // ~>410kms/15A for any R
  - Note: Good for metallic lines with any rotation //
          Very good for H/He lines with any rotation
  - *Note: 'vrg_Z' should only be used if the line is complex and somehow
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
