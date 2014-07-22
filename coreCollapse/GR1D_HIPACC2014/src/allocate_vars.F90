! -*-f90-*-
subroutine allocate_vars

  use GR1D_module

  allocate(x1(n1))
  allocate(x1i(n1))

  allocate(rho(n1))
  allocate(rhop(n1))
  allocate(rhom(n1))

  allocate(v1(n1))
  allocate(v1p(n1))
  allocate(v1m(n1))

  if(do_rotation) then
     allocate(vphi1(n1))
     allocate(vphi1p(n1))
     allocate(vphi1m(n1))
     allocate(vphi(n1))
     allocate(vphip(n1))
     allocate(vphim(n1))
     allocate(omega(n1))
  endif
  allocate(ToverW(n1))

  allocate(v(n1))
  allocate(vp(n1))
  allocate(vm(n1))

  allocate(eps(n1))
  allocate(epsp(n1))
  allocate(epsm(n1))

  allocate(eps_kin(n1))
  allocate(binding_energy(n1))

  allocate(mass(n1))
  allocate(mass1(n1))
  allocate(volume(n1))

  allocate(press(n1))
  allocate(pressth(n1))
  allocate(pressp(n1))
  allocate(pressm(n1))
  allocate(press_nu(n1))
  allocate(dnupdr(n1))
  
  allocate(cs2(n1))
  allocate(cs2p(n1))
  allocate(cs2m(n1))
  
  allocate(temp(n1))

  allocate(entropy(n1))
  allocate(nuchem(n1))
  allocate(massfrac_p(n1))
  allocate(massfrac_n(n1))
  allocate(massfrac_h(n1))
  allocate(massfrac_a(n1))
  allocate(massfrac_abar(n1))
  allocate(massfrac_zbar(n1))

  allocate(ye(n1))
  allocate(yep(n1))
  allocate(yem(n1))

  allocate(q(n1,n_cons))
  allocate(qold(n1,n_cons))
  allocate(qp(n1,n_cons))
  allocate(qm(n1,n_cons))

  allocate(q_hat(n1,n_cons))
  allocate(q_hat_old(n1,n_cons))

  allocate(sqrt_gamma(n1))

  allocate(flux_diff(n1,n_cons))
  allocate(gravsource(n1,n_cons))
  allocate(presssource(n1,n_cons))
  allocate(coolingsource(n1,n_cons))
  allocate(denergyloss(n1))
  
  allocate(atmo(n1))
  
! #############################################
! GR VARIABLES

  allocate(alp(n1))
  allocate(alpp(n1))
  allocate(alpm(n1))
  allocate(phi(n1))
  allocate(phii(n1))
  allocate(dphidr(n1))
  allocate(X(n1))
  allocate(Xp(n1))
  allocate(Xm(n1))
  allocate(W(n1))
  allocate(Wp(n1))
  allocate(Wm(n1))
  allocate(mgrav(n1))
  allocate(mgravi(n1))
  
end subroutine allocate_vars

subroutine deallocate_vars

  use GR1D_module
  implicit none

  deallocate(x1)
  deallocate(x1i)

  deallocate(rho)
  deallocate(rhop)
  deallocate(rhom)

  deallocate(v1)
  deallocate(v1p)
  deallocate(v1m)

  if(do_rotation) then
     deallocate(vphi1)
     deallocate(vphi1p)
     deallocate(vphi1m)
     deallocate(vphi)
     deallocate(vphip)
     deallocate(vphim)
     deallocate(omega)
  endif
  deallocate(ToverW)

  deallocate(v)
  deallocate(vp)
  deallocate(vm)

  deallocate(eps)
  deallocate(epsp)
  deallocate(epsm)

  deallocate(eps_kin)
  deallocate(binding_energy)

  deallocate(mass)
  deallocate(mass1)
  deallocate(volume)

  deallocate(press)
  deallocate(pressth)
  deallocate(pressp)
  deallocate(pressm)
  deallocate(press_nu)
  deallocate(dnupdr)
  
  deallocate(cs2)
  deallocate(cs2p)
  deallocate(cs2m)
  
  deallocate(temp)

  deallocate(entropy)
  deallocate(nuchem)
  deallocate(massfrac_p)
  deallocate(massfrac_n)
  deallocate(massfrac_h)
  deallocate(massfrac_a)
  deallocate(massfrac_abar)
  deallocate(massfrac_zbar)

  deallocate(ye)
  deallocate(yep)
  deallocate(yem)

  deallocate(q)
  deallocate(qold)
  deallocate(qp)
  deallocate(qm)

  deallocate(q_hat)
  deallocate(q_hat_old)

  deallocate(sqrt_gamma)

  deallocate(flux_diff)
  deallocate(gravsource)
  deallocate(presssource)
  deallocate(coolingsource)
  deallocate(denergyloss)
  
  deallocate(atmo)
  
! GR VARIABLES

  deallocate(alp)
  deallocate(alpp)
  deallocate(alpm)
  deallocate(phi)
  deallocate(phii)
  deallocate(dphidr)
  deallocate(X)
  deallocate(Xp)
  deallocate(Xm)
  deallocate(W)
  deallocate(Wp)
  deallocate(Wm)
  deallocate(mgrav)
  deallocate(mgravi)

end subroutine deallocate_vars
