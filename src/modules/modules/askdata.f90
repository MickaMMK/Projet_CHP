module askdatamod

  implicit none

contains

  subroutine askdata(ask)

    implicit none

    integer, intent(out) :: ask
    character(len=1) :: charask

    call getarg(1, charask)
    if(charask == "") then
       print*, "================================================================"
       print*, "------------------ ERROR - Argument manquant -------------------"
       print*, "================================================================"
       print*, "------------ './run 0' pour lire le fichier d'input ------------"
       print*, "- './run 1' pour donner les paramètres au début de l'exécution -"
       print*, "================================================================"
       stop
    else
       read(charask,*) ask
       if(ask < 0 .OR. ask > 1) then
          print*, "================================================================"
          print*, "------------------ ERROR - Argument invalide -------------------"
          print*, "================================================================"
          print*, "------------ './run 0' pour lire le fichier d'input ------------"
          print*, "- './run 1' pour donner les paramètres au début de l'exécution -"
          print*, "================================================================"
          stop
       end if
    end if

  end subroutine askdata

end module askdatamod
