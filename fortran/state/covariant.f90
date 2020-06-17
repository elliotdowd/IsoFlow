
!-*- f90 -*- -
subroutine calc_covariant(s_proj, u, v, U, V)
! calculate covariant velocity for a given cell
! input projected cell face areas [S_zeta.x, S_zeta.y, S_eta.x, S_eta.y, S_zeta, S_eta] and velocity components

    U(:,:) = (1./s_proj(:,:,5)) .* &
           &    (u.*s_proj(:,:,1) + v.*cells.projFaceArea(:,:,2));
       
    V(:,:) = (1./s_proj(:,:,6)) .* &
           &    (u.*s_proj(:,:,3) + v.*s_proj(:,:,4));

end