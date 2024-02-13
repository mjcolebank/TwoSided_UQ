subroutine impedance_driver(tmstps,Period,rho,mu,r_root,r_min,y11,y12,y21,y22,Lr,q,g,&
                            fa1,fa2,fa3,fv1,fv2,fv3,alpha_b,beta_b,lrrA,lrrV,term_ID)
  use f90_tools
  use new_match 
  implicit none




  integer, intent(in)    :: tmstps, term_ID
  real(lng), intent(in)  :: Period,rho,mu,r_root,r_min,Lr,q,g,fa1,fa2,fa3,&
                            fv1,fv2,fv3,alpha_b,beta_b,lrrA,lrrV
  real(lng), intent(out) :: y11(tmstps),y12(tmstps),y21(tmstps),y22(tmstps)
  integer :: j


  do j = 1, tmstps
    y11(j) = 0.0
    y12(j) = 0.0
    y21(j) = 0.0
    y22(j) = 0.0
  end do
  call impedance (tmstps,Period,rho,mu,r_root,r_min,y11,y12,y21,y22,Lr,q,g,&
                    fa1,fa2,fa3,fv1,fv2,fv3,alpha_b,beta_b,lrrA,lrrV,term_ID)



end subroutine impedance_driver
