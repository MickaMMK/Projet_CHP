module getdatamod

  implicit none

contains

  subroutine getdata(N, meth, reproj, pos, transi, raff, remaill, nbp, npart_uni, npart, raff_num, raff_size)

    implicit none

    integer, intent(in) :: N
    integer, intent(out) :: meth, reproj, pos, transi, raff, remaill, nbp, npart, raff_num
    integer, dimension(:), allocatable, intent(out) :: npart_uni
    real(8), intent(out) :: raff_size

    nbp = 0

    call askwarning_int(meth, "Choix de la méthode", (/arg("Eulerien"), arg("Lagrangien")/), 1, 2,&
         & "Erreur dans le choix de la méthode", 1, 2, "", .false.)
    if(meth == 2) then
       call askwarning_int(reproj, "Souhaitez-vous reprojeter les points lagrangiens sur la grille eulerienne ?",&
            & (/arg("Oui"), arg("Non")/), 1, 2, "Erreur dans le choix de la projection", 1, 2, "", .false.)
       if(reproj == 2) then
          reproj = -1
       else if(reproj == 1) then
          call askwarning_int(reproj, "Combien d'itérations entre chaque reprojection ?", (/""/), 1, 1000,&
               & "Erreur dans le choix de la projection", 20, 100, "La valeur donnée pour le nombre d'itération&
               & semble", .true.)
       else
          call abort()
       end if
    end if
    call askwarning_int(pos, "Souhaitez-vous mettre l'eau :", (/arg("Au centre"), arg("A l'extérieur")/),&
         & 1, 2, "Erreur dans le choix de la position de l'eau", 1, 2, "", .false.)
    if(pos == 1) then
       pos = -1
    elseif(pos == 2) then
       pos = 1
    else
       call abort()
    end if
    call askwarning_int(transi, "Quelle transition à l'interface ?", (/arg("Discontinuité"), arg("Linéaire"),&
         & arg("Exponentielle")/), 1, 3, "Erreur dans le choix de la transition à l'interface", 1, 3, "", .false.)
    call askwarning_int(raff, "Souhaitez-vous rajouter des points lagrangiens aux abords de l'interface :",&
         & (/arg("Oui"), arg("Non")/), 1, 2, "Erreur dans le choix du raffinement", 1, 2, "", .false.)
    if(meth == 2) then
       nbp = (N-1)*(N-1)
    end if

    remaill = 2

    if(raff == 1) then
       call askwarning_int(raff_num, "Choisissez le nombre d'anneaux de chaque côté de l'interface &
            &(0 pour du simple front tracking) :", (/""/), 0, 20, "Erreur dans le nombre d'anneaux", 0, 3,&
            & "La valeur donnée pour le nombre d'anneaux semble", .true.)
       allocate(npart_uni(raff_num*2+1))
       call askwarning_int(npart_uni(1), "Choisissez le nombre de points lagrangiens par anneau :", (/""/),&
            & 1, 10000, "Erreur dans le choix du nombre de points", 10, 1000, "La valeur donnée pour le nombre &
            &de points semble", .true.)
       npart_uni = npart_uni(1)
       npart = npart_uni(1)*(raff_num*2+1)
       nbp = nbp + npart
       call askwarning_real(raff_size, "Choisissez la distance entre chaque anneau :", 0.d0, 1.d0,&
            & "Erreur dans le choix de la distance entre chaques anneaux", 1.d-3, 5.d-2, "La valeur donnée pour la&
            & distance semble")
       call askwarning_int(remaill, "Voulez-vous remailler entre chaque itération ?", (/arg("Oui"), arg("Non")/),&
            & 1, 2, "Erreur dans le choix du remaillage", 1, 2, "", .false.)
    end if

  end subroutine getdata

  subroutine askwarning_int(data, ask_sentence, choices_array, min, max, error_sentence, minc, maxc, warning_sentence, warning)

    implicit none

    integer, intent(out) :: data
    integer, intent(in) :: min, max, minc, maxc
    character(len=*), intent(in) :: ask_sentence, error_sentence, warning_sentence
    character(len=*), dimension(:), intent(in) :: choices_array
    logical, intent(in) :: warning
    logical :: ok
    character(len=100) :: charmin, charmax, charminc, charmaxc

    write(charmin,*) min
    write(charmax,*) max
    write(charminc,*) minc
    write(charmaxc,*) maxc


    ok = .false.
    do while(ok .EQV. .false.)
       call askdata_int(data, ask_sentence, choices_array)
       if(warning) then
          call checkdata_int(data, min, max, error_sentence//" -- min : "//trim(adjustl(charmin))//", max : "&
               & //trim(adjustl(charmax))//" -- Range proposée : "//trim(adjustl(charminc))//" -> "&
               & //trim(adjustl(charmaxc)), ok)
          if(ok .EQV. .true.) then
             call warning_int(data, minc, max, warning_sentence//" faible -- Range proposée : "&
                  & //trim(adjustl(charminc))//" -> "//trim(adjustl(charmaxc)), ok)
          end if
          if(ok .EQV. .true.) then
             call warning_int(data, min, maxc, warning_sentence//" élevée -- Range proposée : "&
                  & //trim(adjustl(charminc))//" -> "//trim(adjustl(charmaxc)), ok)
          end if
       else
          call checkdata_int(data, min, max, error_sentence//" -- min : "//trim(adjustl(charmin))//", max : "&
               & //trim(adjustl(charmax)), ok)
       end if
    end do

  end subroutine askwarning_int

  subroutine askwarning_real(data, ask_sentence, min, max, error_sentence, minc, maxc, warning_sentence)

    implicit none

    real(8), intent(out) :: data
    real(8), intent(in) :: min, max, minc, maxc
    character(len=*), intent(in) :: ask_sentence, error_sentence, warning_sentence
    logical :: ok
    character(len=100) :: charmin, charmax, charminc, charmaxc

    write(charmin,*) min
    write(charmax,*) max
    write(charminc,*) minc
    write(charmaxc,*) maxc


    ok = .false.
    do while(ok .EQV. .false.)
       call askdata_real(data, ask_sentence)
       call checkdata_real(data, min, max, error_sentence//" -- min : "//trim(adjustl(charmin))//", max : "&
            & //trim(adjustl(charmax))//" -- Range proposée : "//trim(adjustl(charminc))//" -> "&
            & //trim(adjustl(charmaxc)), ok)
       if(ok .EQV. .true.) then
          call warning_real(data, minc, max, warning_sentence//" faible -- Range proposée : "&
               & //trim(adjustl(charminc))//" -> "//trim(adjustl(charmaxc)), ok)
       end if
       if(ok .EQV. .true.) then
          call warning_real(data, min, maxc, warning_sentence//" élevée -- Range proposée : "&
               & //trim(adjustl(charminc))//" -> "//trim(adjustl(charmaxc)), ok)
       end if
    end do

  end subroutine askwarning_real

  subroutine askdata_int(data, ask_sentence, choices_array)

    implicit none

    integer, intent(out) :: data
    character(len=*), intent(in) :: ask_sentence
    character(len=*), dimension(:) :: choices_array
    integer :: i
    character(len=1) :: num

    print*, ask_sentence
    if(size(choices_array) > 1) then
       do i = 1, size(choices_array)
          write(num, '(I1)') i
          print*, num//" - "//trim(choices_array(i))
       end do
    end if
    read*, data

  end subroutine askdata_int

  subroutine askdata_real(data, ask_sentence)

    implicit none

    real(8), intent(out) :: data
    character(len=*), intent(in) :: ask_sentence

    print*, ask_sentence
    read*, data

  end subroutine askdata_real

  subroutine checkdata_int(data, min, max, error_sentence, ok)

    implicit none

    integer, intent(in) :: data, min, max
    character(len=*), intent(in) :: error_sentence
    logical, intent(inout) :: ok
    
    if(data < min .OR. data > max) then
       print*, "- ERROR - "//error_sentence
       ok = .false.
    else
       ok = .true.
    end if

  end subroutine checkdata_int

  subroutine checkdata_real(data, min, max, error_sentence, ok)

    implicit none

    real(8), intent(in) :: data, min, max
    character(len=*), intent(in) :: error_sentence
    logical, intent(inout) :: ok

    if(data < min .OR. data > max) then
       print*, "- ERROR - "//error_sentence
       ok = .false.
    else
       ok = .true.
    end if

  end subroutine checkdata_real

  subroutine warning_int(data, min, max, warning_sentence, ok)

    implicit none

    integer, intent(in) :: data, min, max
    character(len=*), intent(in) :: warning_sentence
    logical, intent(inout) :: ok
    integer :: avis
    
    if(data < min .OR. data > max) then
       print*, "- WARNING - "//warning_sentence
       call askwarning_int(avis, "Etes-vous sûr de votre choix ?", (/arg("Continuer"), arg("Changer la valeur")/),&
            & 1, 2, "Erreur dans le choix après warning", 1, 2, "", .false.)
       if(avis == 1) then
          ok = .true.
       elseif(avis == 2) then
          ok = .false.
       else
          call abort()
       end if
    else
       ok = .true.
    end if

  end subroutine warning_int

  subroutine warning_real(data, min, max, warning_sentence, ok)

    implicit none

    real(8), intent(in) :: data, min, max
    character(len=*), intent(in) :: warning_sentence
    logical, intent(inout) :: ok
    integer :: avis
    
    if(data < min .OR. data > max) then
       print*, "- WARNING - "//warning_sentence
       call askdata_int(avis, "Etes-vous sûr de votre choix ?", (/arg("Continuer"), arg("Changer la valeur")/))
       call checkdata_int(avis, 1, 2, "Erreur dans le choix après warning", ok)
       if(avis == 1) then
          ok = .true.
       elseif(avis == 2) then
          ok = .false.
       else
          call abort()
       end if
    else
       ok = .true.
    end if

  end subroutine warning_real

  function arg(word) result(big_word)

    implicit none

    character(len=*), intent(in) :: word
    character(len=100) :: big_word

    big_word = word

  end function arg

  subroutine abort()

    implicit none

    print*, "=== UNKNOWN ERROR IN GETDATA.F90 ==="
    print*, "===================================="
    print*, "============== ABORT ==============="
    print*, "===================================="
    stop

  end subroutine abort

end module getdatamod
