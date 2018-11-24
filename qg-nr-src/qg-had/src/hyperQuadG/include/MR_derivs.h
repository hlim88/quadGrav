        call deriv_x(q(dx_alpha),alpha,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_alpha),alpha,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_alpha),alpha,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'alpha'
          call check_derivs_1(rc, alpha, q(dx_alpha),
     &                      q(dy_alpha),q(dz_alpha),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_shiftx),shiftx,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_shiftx),shiftx,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_shiftx),shiftx,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'shiftx'
          call check_derivs_1(rc, shiftx, q(dx_shiftx),
     &                      q(dy_shiftx),q(dz_shiftx),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_shifty),shifty,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_shifty),shifty,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_shifty),shifty,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'shifty'
          call check_derivs_1(rc, shifty, q(dx_shifty),
     &                      q(dy_shifty),q(dz_shifty),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_shiftz),shiftz,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_shiftz),shiftz,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_shiftz),shiftz,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'shiftz'
          call check_derivs_1(rc, shiftz, q(dx_shiftz),
     &                      q(dy_shiftz),q(dz_shiftz),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_gbx),gbx,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_gbx),gbx,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_gbx),gbx,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gbx'
          call check_derivs_1(rc, gbx, q(dx_gbx),
     &                      q(dy_gbx),q(dz_gbx),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_gby),gby,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_gby),gby,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_gby),gby,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gby'
          call check_derivs_1(rc, gby, q(dx_gby),
     &                      q(dy_gby),q(dz_gby),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_gbz),gbz,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_gbz),gbz,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_gbz),gbz,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gbz'
          call check_derivs_1(rc, gbz, q(dx_gbz),
     &                      q(dy_gbz),q(dz_gbz),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_chi),chi,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_chi),chi,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_chi),chi,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'chi'
          call check_derivs_1(rc, chi, q(dx_chi),
     &                      q(dy_chi),q(dz_chi),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Gamtx),Gamtx,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Gamtx),Gamtx,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Gamtx),Gamtx,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Gamtx'
          call check_derivs_1(rc, Gamtx, q(dx_Gamtx),
     &                      q(dy_Gamtx),q(dz_Gamtx),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Gamty),Gamty,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Gamty),Gamty,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Gamty),Gamty,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Gamty'
          call check_derivs_1(rc, Gamty, q(dx_Gamty),
     &                      q(dy_Gamty),q(dz_Gamty),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Gamtz),Gamtz,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Gamtz),Gamtz,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Gamtz),Gamtz,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Gamtz'
          call check_derivs_1(rc, Gamtz, q(dx_Gamtz),
     &                      q(dy_Gamtz),q(dz_Gamtz),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_trK),trK,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_trK),trK,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_trK),trK,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'trK'
          call check_derivs_1(rc, trK, q(dx_trK),
     &                      q(dy_trK),q(dz_trK),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_gtxx),gtxx,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtxx),gtxx,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtxx),gtxx,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtxx'
          call check_derivs_1(rc, gtxx, q(dx_gtxx),
     &                      q(dy_gtxx),q(dz_gtxx),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_gtxy),gtxy,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtxy),gtxy,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtxy),gtxy,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtxy'
          call check_derivs_1(rc, gtxy, q(dx_gtxy),
     &                      q(dy_gtxy),q(dz_gtxy),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_gtxz),gtxz,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtxz),gtxz,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtxz),gtxz,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtxz'
          call check_derivs_1(rc, gtxz, q(dx_gtxz),
     &                      q(dy_gtxz),q(dz_gtxz),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_gtyy),gtyy,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtyy),gtyy,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtyy),gtyy,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtyy'
          call check_derivs_1(rc, gtyy, q(dx_gtyy),
     &                      q(dy_gtyy),q(dz_gtyy),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_gtyz),gtyz,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtyz),gtyz,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtyz),gtyz,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtyz'
          call check_derivs_1(rc, gtyz, q(dx_gtyz),
     &                      q(dy_gtyz),q(dz_gtyz),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_gtzz),gtzz,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_gtzz),gtzz,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_gtzz),gtzz,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtzz'
          call check_derivs_1(rc, gtzz, q(dx_gtzz),
     &                      q(dy_gtzz),q(dz_gtzz),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Atxx),Atxx,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Atxx),Atxx,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Atxx),Atxx,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Atxx'
          call check_derivs_1(rc, Atxx, q(dx_Atxx),
     &                      q(dy_Atxx),q(dz_Atxx),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Atxy),Atxy,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Atxy),Atxy,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Atxy),Atxy,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Atxy'
          call check_derivs_1(rc, Atxy, q(dx_Atxy),
     &                      q(dy_Atxy),q(dz_Atxy),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Atxz),Atxz,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Atxz),Atxz,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Atxz),Atxz,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Atxz'
          call check_derivs_1(rc, Atxz, q(dx_Atxz),
     &                      q(dy_Atxz),q(dz_Atxz),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Atyy),Atyy,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Atyy),Atyy,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Atyy),Atyy,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Atyy'
          call check_derivs_1(rc, Atyy, q(dx_Atyy),
     &                      q(dy_Atyy),q(dz_Atyy),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Atyz),Atyz,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Atyz),Atyz,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Atyz),Atyz,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Atyz'
          call check_derivs_1(rc, Atyz, q(dx_Atyz),
     &                      q(dy_Atyz),q(dz_Atyz),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Atzz),Atzz,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Atzz),Atzz,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Atzz),Atzz,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Atzz'
          call check_derivs_1(rc, Atzz, q(dx_Atzz),
     &                      q(dy_Atzz),q(dz_Atzz),
     &                   fname, deriv_time, myproc,
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Ex),Ex,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Ex),Ex,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Ex),Ex,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Ex'
          call check_derivs_1(rc, Ex, q(dx_Ex),
     &                    q(dy_Ex),q(dz_Ex),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Ey),Ey,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Ey),Ey,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Ey),Ey,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Ey'
          call check_derivs_1(rc, Ey, q(dx_Ey),
     &                    q(dy_Ey),q(dz_Ey),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Ez),Ez,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Ez),Ez,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Ez),Ez,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Ez'
          call check_derivs_1(rc, Ez, q(dx_Ez),
     &                    q(dy_Ez),q(dz_Ez),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Bx),Bx,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Bx),Bx,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Bx),Bx,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Bx'
          call check_derivs_1(rc, Bx, q(dx_Bx),
     &                    q(dy_Bx),q(dz_Bx),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_By),By,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_By),By,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_By),By,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'By'
          call check_derivs_1(rc, By, q(dx_By),
     &                    q(dy_By),q(dz_By),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Bz),Bz,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Bz),Bz,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Bz),Bz,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Bz'
          call check_derivs_1(rc, Bz, q(dx_Bz),
     &                    q(dy_Bz),q(dz_Bz),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Phi_em),Phi_em,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Phi_em),Phi_em,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Phi_em),Phi_em,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Phi_em'
          call check_derivs_1(rc, Phi_em, q(dx_Phi_em),
     &                    q(dy_Phi_em),q(dz_Phi_em),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_Psi_em),Psi_em,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_Psi_em),Psi_em,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_Psi_em),Psi_em,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'Psi_em'
          call check_derivs_1(rc, Psi_em, q(dx_Psi_em),
     &                    q(dy_Psi_em),q(dz_Psi_em),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_phiR),phiR,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_phiR),phiR,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_phiR),phiR,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'phiR'
          call check_derivs_1(rc, phiR, q(dx_phiR),
     &                    q(dy_phiR),q(dz_phiR),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_phiI),phiI,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_phiI),phiI,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_phiI),phiI,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'phiI'
          call check_derivs_1(rc, phiI, q(dx_phiI),
     &                    q(dy_phiI),q(dz_phiI),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_piR),piR,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_piR),piR,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_piR),piR,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'piR'
          call check_derivs_1(rc, piR, q(dx_piR),
     &                    q(dy_piR),q(dz_piR),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_x(q(dx_piI),piI,dx, nx, ny, nz, d_type)
        call deriv_y(q(dy_piI),piI,dy, nx, ny, nz, d_type)
        call deriv_z(q(dz_piI),piI,dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'piI'
          call check_derivs_1(rc, piI, q(dx_piI),
     &                    q(dy_piI),q(dz_piI),
     &                 fname, deriv_time, myproc,
     &                 x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_gtxx),gtxx,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_gtxx),gtxx,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_gtxx),gtxx,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_gtxx),q(dx_gtxx),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_gtxx),q(dx_gtxx),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_gtxx),q(dy_gtxx),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtxx'
          call check_derivs_2(rc, gtxx,
     &                   q(dxx_gtxx),q(dxy_gtxx),
     &                   q(dxz_gtxx),q(dyy_gtxx),
     &                   q(dyz_gtxx),q(dzz_gtxx),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_gtxy),gtxy,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_gtxy),gtxy,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_gtxy),gtxy,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_gtxy),q(dx_gtxy),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_gtxy),q(dx_gtxy),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_gtxy),q(dy_gtxy),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtxy'
          call check_derivs_2(rc, gtxy,
     &                   q(dxx_gtxy),q(dxy_gtxy),
     &                   q(dxz_gtxy),q(dyy_gtxy),
     &                   q(dyz_gtxy),q(dzz_gtxy),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_gtxz),gtxz,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_gtxz),gtxz,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_gtxz),gtxz,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_gtxz),q(dx_gtxz),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_gtxz),q(dx_gtxz),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_gtxz),q(dy_gtxz),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtxz'
          call check_derivs_2(rc, gtxz,
     &                   q(dxx_gtxz),q(dxy_gtxz),
     &                   q(dxz_gtxz),q(dyy_gtxz),
     &                   q(dyz_gtxz),q(dzz_gtxz),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_gtyy),gtyy,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_gtyy),gtyy,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_gtyy),gtyy,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_gtyy),q(dx_gtyy),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_gtyy),q(dx_gtyy),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_gtyy),q(dy_gtyy),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtyy'
          call check_derivs_2(rc, gtyy,
     &                   q(dxx_gtyy),q(dxy_gtyy),
     &                   q(dxz_gtyy),q(dyy_gtyy),
     &                   q(dyz_gtyy),q(dzz_gtyy),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_gtyz),gtyz,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_gtyz),gtyz,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_gtyz),gtyz,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_gtyz),q(dx_gtyz),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_gtyz),q(dx_gtyz),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_gtyz),q(dy_gtyz),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtyz'
          call check_derivs_2(rc, gtyz,
     &                   q(dxx_gtyz),q(dxy_gtyz),
     &                   q(dxz_gtyz),q(dyy_gtyz),
     &                   q(dyz_gtyz),q(dzz_gtyz),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_gtzz),gtzz,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_gtzz),gtzz,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_gtzz),gtzz,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_gtzz),q(dx_gtzz),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_gtzz),q(dx_gtzz),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_gtzz),q(dy_gtzz),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'gtzz'
          call check_derivs_2(rc, gtzz,
     &                   q(dxx_gtzz),q(dxy_gtzz),
     &                   q(dxz_gtzz),q(dyy_gtzz),
     &                   q(dyz_gtzz),q(dzz_gtzz),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_chi),chi,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_chi),chi,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_chi),chi,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_chi),q(dx_chi),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_chi),q(dx_chi),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_chi),q(dy_chi),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'chi'
          call check_derivs_2(rc, chi,
     &                   q(dxx_chi),q(dxy_chi),
     &                   q(dxz_chi),q(dyy_chi),
     &                   q(dyz_chi),q(dzz_chi),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_alpha),alpha,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_alpha),alpha,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_alpha),alpha,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_alpha),q(dx_alpha),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_alpha),q(dx_alpha),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_alpha),q(dy_alpha),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'alpha'
          call check_derivs_2(rc, alpha,
     &                   q(dxx_alpha),q(dxy_alpha),
     &                   q(dxz_alpha),q(dyy_alpha),
     &                   q(dyz_alpha),q(dzz_alpha),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_shiftx),shiftx,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_shiftx),shiftx,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_shiftx),shiftx,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_shiftx),q(dx_shiftx),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_shiftx),q(dx_shiftx),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_shiftx),q(dy_shiftx),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'shiftx'
          call check_derivs_2(rc, shiftx,
     &                   q(dxx_shiftx),q(dxy_shiftx),
     &                   q(dxz_shiftx),q(dyy_shiftx),
     &                   q(dyz_shiftx),q(dzz_shiftx),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_shifty),shifty,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_shifty),shifty,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_shifty),shifty,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_shifty),q(dx_shifty),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_shifty),q(dx_shifty),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_shifty),q(dy_shifty),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'shifty'
          call check_derivs_2(rc, shifty,
     &                   q(dxx_shifty),q(dxy_shifty),
     &                   q(dxz_shifty),q(dyy_shifty),
     &                   q(dyz_shifty),q(dzz_shifty),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_shiftz),shiftz,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_shiftz),shiftz,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_shiftz),shiftz,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_shiftz),q(dx_shiftz),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_shiftz),q(dx_shiftz),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_shiftz),q(dy_shiftz),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'shiftz'
          call check_derivs_2(rc, shiftz,
     &                   q(dxx_shiftz),q(dxy_shiftz),
     &                   q(dxz_shiftz),q(dyy_shiftz),
     &                   q(dyz_shiftz),q(dzz_shiftz),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_phiR),phiR,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_phiR),phiR,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_phiR),phiR,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_phiR),q(dx_phiR),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_phiR),q(dx_phiR),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_phiR),q(dy_phiR),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'phiR'
          call check_derivs_2(rc, phiR,
     &                   q(dxx_phiR),q(dxy_phiR),
     &                   q(dxz_phiR),q(dyy_phiR),
     &                   q(dyz_phiR),q(dzz_phiR),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
        call deriv_xx(q(dxx_phiI),phiI,dx, nx, ny, nz, dd_type)
        call deriv_yy(q(dyy_phiI),phiI,dy, nx, ny, nz, dd_type)
        call deriv_zz(q(dzz_phiI),phiI,dz, nx, ny, nz, dd_type)
        call deriv_y(q(dxy_phiI),q(dx_phiI),dy, nx, ny, nz, d_type)
        call deriv_z(q(dxz_phiI),q(dx_phiI),dz, nx, ny, nz, d_type)
        call deriv_z(q(dyz_phiI),q(dy_phiI),dz, nx, ny, nz, d_type)
        if (debug_derivs) then
          fname = 'phiI'
          call check_derivs_2(rc, phiI,
     &                   q(dxx_phiI),q(dxy_phiI),
     &                   q(dxz_phiI),q(dyy_phiI),
     &                   q(dyz_phiI),q(dzz_phiI),
     &                   fname, deriv_time, myproc, 
     &                   x1d, y1d, z1d, nx, ny, nz)
          if (rc .eq. 0) deriv_error = 1
        end if
