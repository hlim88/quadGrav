        call adv_deriv_x(q(adx_gtxx),gtxx,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_gtxx),gtxx,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_gtxx),gtxx,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_gtxy),gtxy,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_gtxy),gtxy,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_gtxy),gtxy,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_gtxz),gtxz,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_gtxz),gtxz,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_gtxz),gtxz,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_gtyy),gtyy,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_gtyy),gtyy,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_gtyy),gtyy,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_gtyz),gtyz,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_gtyz),gtyz,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_gtyz),gtyz,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_gtzz),gtzz,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_gtzz),gtzz,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_gtzz),gtzz,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_Atxx),Atxx,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_Atxx),Atxx,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_Atxx),Atxx,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_Atxy),Atxy,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_Atxy),Atxy,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_Atxy),Atxy,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_Atxz),Atxz,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_Atxz),Atxz,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_Atxz),Atxz,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_Atyy),Atyy,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_Atyy),Atyy,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_Atyy),Atyy,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_Atyz),Atyz,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_Atyz),Atyz,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_Atyz),Atyz,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_Atzz),Atzz,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_Atzz),Atzz,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_Atzz),Atzz,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_alpha),alpha,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_alpha),alpha,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_alpha),alpha,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_shiftx),shiftx,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_shiftx),shiftx,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_shiftx),shiftx,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_shifty),shifty,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_shifty),shifty,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_shifty),shifty,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_shiftz),shiftz,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_shiftz),shiftz,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_shiftz),shiftz,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_chi),chi,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_chi),chi,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_chi),chi,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_Gamtx),Gamtx,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_Gamtx),Gamtx,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_Gamtx),Gamtx,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_Gamty),Gamty,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_Gamty),Gamty,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_Gamty),Gamty,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_Gamtz),Gamtz,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_Gamtz),Gamtz,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_Gamtz),Gamtz,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_trK),trK,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_trK),trK,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_trK),trK,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_gbx),gbx,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_gbx),gbx,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_gbx),gbx,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_gby),gby,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_gby),gby,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_gby),gby,dz, nx, ny, nz,
     &          d_type, shiftz)
        call adv_deriv_x(q(adx_gbz),gbz,dx, nx, ny, nz,
     &          d_type, shiftx)
        call adv_deriv_y(q(ady_gbz),gbz,dy, nx, ny, nz,
     &          d_type, shifty)
        call adv_deriv_z(q(adz_gbz),gbz,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_Ex),Ex,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_Ex),Ex,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_Ex),Ex,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_Ey),Ey,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_Ey),Ey,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_Ey),Ey,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_Ez),Ez,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_Ez),Ez,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_Ez),Ez,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_Bx),Bx,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_Bx),Bx,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_Bx),Bx,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_By),By,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_By),By,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_By),By,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_Bz),Bz,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_Bz),Bz,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_Bz),Bz,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_Phi_em),Phi_em,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_Phi_em),Phi_em,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_Phi_em),Phi_em,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_Psi_em),Psi_em,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_Psi_em),Psi_em,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_Psi_em),Psi_em,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_phiR),phiR,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_phiR),phiR,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_phiR),phiR,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_phiI),phiI,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_phiI),phiI,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_phiI),phiI,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_piR),piR,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_piR),piR,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_piR),piR,dz, nx, ny, nz,
     &          d_type, shiftz)
          call adv_deriv_x(q(adx_piI),piI,dx, nx, ny, nz,
     &          d_type, shiftx)
          call adv_deriv_y(q(ady_piI),piI,dy, nx, ny, nz,
     &          d_type, shifty)
          call adv_deriv_z(q(adz_piI),piI,dz, nx, ny, nz,
     &          d_type, shiftz)
