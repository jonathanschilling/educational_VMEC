!> \file
MODULE realspace

  USE stel_kinds

  IMPLICIT NONE

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: r1 !< \f$R\f$
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ru !< \f$\partial R/\partial \theta\f$
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rv !< \f$\partial R/\partial \zeta\f$
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: z1 !< \f$Z \f$
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: zu   !< \f$\partial Z/\partial \theta\f$
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: zv   !< \f$\partial Z/\partial \zeta\f$

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rcon !< spectral condensation term in \f$R\f$
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: zcon !< spectral condensation term in \f$Z\f$
  REAL(rprec), DIMENSION(:),   ALLOCATABLE :: rcon0 !< spectral condensation term in \f$R\f$ at start of current multi-grid iteration
  REAL(rprec), DIMENSION(:),   ALLOCATABLE :: zcon0 !< spectral condensation term in \f$R\f$ at start of current multi-grid iteration

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: guu !< \f$g_{\theta, \theta}\f$ metric element
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: guv !< \f$g_{\theta, \zeta}\f$ metric element
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: gvv !< \f$g_{\zeta, \zeta}\f$ metric element

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: ru0 !< \f$\partial R/\partial \theta\f$, even-m and odd-m added together appropriately
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: zu0 !< \f$\partial Z/\partial \theta\f$, even-m and odd-m added together appropriately

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: gcon !< spectral condensation force; "alias force"

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: phip  !< radial derivative of phi/(2*pi) on half-grid
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: chip  !< radial derivative of chi/(2*pi) on half-grid
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: shalf !< sqrt(s) ,two-dimensional array on half-grid
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: sqrts !< sqrt(s), two-dimensional array on full-grid

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: wint  !< two-dimensional array for normalizing angle integrations

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: extra1
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: extra2
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: extra3
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: extra4

END MODULE realspace
