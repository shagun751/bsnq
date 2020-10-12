## Vertical velocity angular

1. [Wave-angle [2020-07-03]](#log_bsnqM_v0004_1)
1. [VertVelAng Attempt-1 [2020-07-05]](#log_bsnqM_v0004_2)


### Attempting
- Calculate velocities along the depth for angular waves
- Try and avoid calculation of 'w' velocity. We already avoided need for pressure


### List of Work
- [x] Calculate local wave-angle. _locWvAng()_
- [x] BUG_FIX: wavHReset had a mistake. <br>Only mod(12.01,3) would work. mod(11.99,3) wont work. Corrected.
- [ ] Angular wave expression for depth-resolved vel.

-----------------------------------------------


<a name = 'log_bsnqM_v0004_1' ></a>

### VertVelAng Attempt-1 [2020-07-05]
- Not using wave-driection _(nx, ny)_ because that becomes too complex and unreliable. Wave direction is notoriously difficult to find.
- Assuming that _u_ only depends on _(p, z, h)_ and _v_ only depends on _(q, z, h)_
    - _u_ and _v_ are calculated using the exact same expression as the uni-directional wave, except that for _v_ we replace _p_ with _q_ and x-derivatives with y-derivatives.
- Velocity _w_ depends on _(p, q, z, h)_.
    - Calcualted the expression for unidirectional wave using _(p, z, h)_ and x-derivatives.
    - Add to the above the expressions for unidirectional wave using _(q, z, h)_ and y-derivatives.
    - No cross-derivatives

The expressions are given below, done in the Matlab code 'OtherCodes/Matlab/VertVel_Angular/vertVel_3.m'

```
    u0c = um0 - 0.5*d * ( umd0_xx ) ...
           + d*d/6 * ( um0_xx ) ...
           - z0 * ( umd0_xx ) ...
           - z0*z0/2 * ( um0_xx );

    v0c = vm0 - 0.5*d * ( vmd0_yy ) ...
           + d*d/6 * ( vm0_yy ) ...
           - z0 * ( vmd0_yy ) ...
           - z0*z0/2 * ( vm0_yy );

    w0c = - ( umd0_x ) ...
          - z0 * ( um0_x ) ...
          + z0*d/2 * ( umd0_xxx ) ...
          - z0*d*d/6 * ( um0_xxx ) ...
          + z0*z0/2 * ( umd0_xxx ) ...
          + z0*z0*z0/6 * ( um0_xxx );
          - ( vmd0_y ) ...
          - z0 * ( vm0_y ) ...
          + z0*d/2 * ( vmd0_yyy ) ...
          - z0*d*d/6 * ( vm0_yyy ) ...
          + z0*z0/2 * ( vmd0_yyy ) ...
          + z0*z0*z0/6 * ( vm0_yyy );
```

-----------------------------------------------


<a name = 'log_bsnqM_v0004_1' ></a>

### Wave-angle [2020-07-03]

#### Method-Vel [mostly wrong]
- The convectional expression for wave angle is using wave-elevation. <br> tan( &theta; ) = (d&eta;/dy) / (d&eta;/dx) as given in Sorenson (2004) while describing roller-breaker.
- But instad of calculating the derivatives I am relying on depth-integ vel P Q and hoping the results tan( &theta; ) = Q / P will be alright
- Calculated if the velMag .gt. 1d-10

#### Method-Eta [mostly wrong]
- The convectional expression for wave angle is using wave-elevation. <br> tan( &theta; ) = (d&eta;/dy) / (d&eta;/dx) as given in Sorenson (2004) while describing roller-breaker.
- Calculated using MLS shape function.
- No visible effect on run-time.

#### Results

| |
| :-------------: |
| **Figure :** Wave-Angle calculated using Vel vs calculated using Eta |
| <img width="90%" src="./log0004/fc45a_wvAng.gif"> |
| **Figure :** Wave-Angle calculated using Vel vs calculated using Eta along centre-line of the domain |
| <img width="90%" src="./log0004/fc45a_wvAng_centreLine.gif"> |
| **Figure :** Wave-Angle calculated using Vel vs calculated using Eta along centre-line of the domain at t = 20s|
| <img width="90%" src="./log0004/wvAngUsing_VelBlk_EtaRed_t20s.png"> |

The observations are :

- The local instantaneous wave angle from Vel and Eta have a phase difference of 90deg for the uni-directional wave, before the collision with the wall.
- After the collision with the wall the returning wave super-imposing on the forward wave give similar wave angle as calculated using Vel and Eta.
- Method-Vel is more noisy and sensitive.
- The &Eta; wave elevation is anyways smoother than velocity, and additionally the MLS derivative ensures its more smooth.
- Method-Eta is mentioned in a paper, hence it may be correct, but neither seem very useful for the particle velocity calculations!
- Keeping the Method-Eta and commenting out the Method-Vel.

**For now I have disabled _lovWvAng()_ inside _postInstructs()_ by commenting it because it is not needed by default.**

-----------------------------------------------

## References

1. Madsen, P. A., & Agnon, Y. (2003). Accuracy and convergence of velocity formulations for water waves in the framework of Boussinesq theory. Journal of Fluid Mechanics, 477(477), 285–319. https://doi.org/10.1017/S0022112002003257