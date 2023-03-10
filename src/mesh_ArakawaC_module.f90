MODULE mesh_ArakawaC_module

  ! Routines used in creating the Arakawa C staggered mesh.

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  
  ! Import specific functionality           
  USE data_types_module,               ONLY: type_mesh
  USE mesh_help_functions_module,      ONLY: is_boundary_segment

  IMPLICIT NONE

  CONTAINS
    
  SUBROUTINE make_Ac_mesh( mesh)
    ! Determine the Arakawa C mesh vertices and their neighbour functions.
    ! For a detailed explanation, read the documentation (not kidding, its in there).    
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'make_Ac_mesh'
    INTEGER                                       :: vi, ci, vj, cj, iti, ti, n1, n2, n3, vil, vir, til, tir
    INTEGER,  DIMENSION(:,:), ALLOCATABLE         :: Aci_temp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: VAc_temp
    
    ! Add routine to path
    CALL init_routine( routine_name)
                
    ! Allocate memory for the mapping arrays. Aci is stored in temporary memory on the master process first,
    ! once we know how much is needed, we will allocate shared memory and transfer the data there.
    if (allocated (mesh%iAci)) deallocate(mesh%iAci)
    allocate( mesh%iAci( mesh%nV, mesh%nC_mem))
    
    ! Go through all vertex connections, determine the location of the Ac vertex (on the connection midpoint),
    ! the indices of the left and right Aa vertices, and the neighbour functions on the Ac vertex
    
  
    ! Allocate temporary memory for VAc, Aci and neighbour functions
    ALLOCATE(VAc_temp( mesh%nTri * 3,2))
    ALLOCATE(Aci_temp( mesh%nTri * 3,6))
    
    mesh%nAc   = 0
    mesh%iAci  = 0
    Aci_temp   = 0
    VAc_temp   = 0._dp
    
    ! Go over all Aa vertices
    DO vi = 1, mesh%nV
    
      ! Go over all the vertex connections
      DO ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        
        ! Skip connections that were already considered in the opposite direction
        IF (mesh%iAci( vi,ci) > 0) CYCLE
        
        ! List Ac vertex number in reference array
        mesh%nAc = mesh%nAc + 1
        mesh%iAci( vi,ci) = mesh%nAc
        VAc_temp(mesh%nAc,:) = (mesh%V( vi,:) + mesh%V(vj,:)) / 2._dp
          
        ! Find reverse connection
        DO cj = 1, mesh%nC( vj)
          IF (mesh%C( vj,cj) == vi) THEN
            mesh%iAci( vj,cj) = mesh%nAC
            EXIT
          END IF
        END DO
        
        IF (.NOT. is_boundary_segment( mesh, vi, vj)) THEN
            
          ! Find indices of vl and vr
          vil = 0
          vir = 0
          til = 0
          tir = 0
          DO iti = 1, mesh%niTri( vi)
            ti = mesh%iTri( vi,iti)
            DO n1 = 1, 3
              n2 = n1+1
              IF (n2 == 4) n2 = 1
              n3 = n2+1
              IF (n3 == 4) n3 = 1
              
              IF     (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj) THEN
                til = ti
                vil = mesh%Tri( ti,n3)
              ELSEIF (mesh%Tri( ti,n1) == vj .AND. mesh%Tri( ti,n2) == vi) THEN
                tir = ti
                vir = mesh%Tri( ti,n3)
              END IF
            END DO
          END DO
          
        ELSE
          ! Surface slope on an edge segment Ac vertex equals that of the (single) triangle
            
          ! Find index of vl
          vil = 0
          vir = 0
          til = 0
          tir = 0
          DO iti = 1, mesh%niTri( vi)
            ti = mesh%iTri( vi,iti)
            DO n1 = 1, 3
              n2 = n1+1
              IF (n2 == 4) n2 = 1
              n3 = n2+1
              IF (n3 == 4) n3 = 1
              
              IF     (mesh%Tri( ti,n1) == vi .AND. mesh%Tri( ti,n2) == vj) THEN
                til = ti
                vil = mesh%Tri( ti,n3)
              ELSEIF (mesh%Tri( ti,n1) == vj .AND. mesh%Tri( ti,n2) == vi) THEN
                tir = ti
                vir = mesh%Tri( ti,n3)
              END IF
        
            END DO
          END DO         
        
        END IF

        ! List relevant Aa vertex indices vor Ac vertex in reference array
        Aci_temp( mesh%nAc,:) = [vi, vj, vil, vir, til, tir]
          
      END DO ! DO ci = 1, mesh%nC(vi)
    END DO ! DO vi = 1, mesh%nV
  
    ! Allocate shared memory, move data there
    if (allocated(mesh%VAc)) deallocate(mesh%VAc)
    if (allocated(mesh%Aci)) deallocate(mesh%Aci)
    allocate(mesh%VAc( mesh%nAc, 2))
    allocate(mesh%Aci( mesh%nAc, 6))
    
    mesh%VAc   = VAc_temp(   1:mesh%nAc,:)
    mesh%Aci   = Aci_temp(   1:mesh%nAc,:)
    
    ! Deallocate temporary memory
    DEALLOCATE( VAc_temp)
    DEALLOCATE( Aci_temp)
    
    ! Determine Ac vertex domains
    CALL partition_list( mesh%nAc, par%i, par%n, mesh%ci1, mesh%ci2)
    
    ! Find Ac edge indices
    CALL find_Ac_edge_indices( mesh)  
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 5)
    
  END SUBROUTINE make_Ac_mesh
  SUBROUTINE find_Ac_edge_indices( mesh)
    ! Find the edge indices of the Ac vertices
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'find_Ac_edge_indices'
    INTEGER                                       :: aci, vi, vj
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate shared memory
    if (allocated(mesh%edge_index_Ac)) deallocate(mesh%edge_index_Ac)
    allocate(mesh%edge_index_Ac(1:mesh%nAc))
    
    DO aci = 1, mesh%nAc ! TODO parallel?
  
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      
      IF     ((mesh%edge_index( vi)==8 .OR. mesh%edge_index( vi) == 1 .OR. mesh%edge_index( vi) == 2) .AND. &
              (mesh%edge_index( vj)==8 .OR. mesh%edge_index( vj) == 1 .OR. mesh%edge_index( vj) == 2)) THEN
        ! North
        mesh%edge_index_Ac( aci) = 1
      ELSEIF ((mesh%edge_index( vi)==2 .OR. mesh%edge_index( vi) == 3 .OR. mesh%edge_index( vi) == 4) .AND. &
              (mesh%edge_index( vj)==2 .OR. mesh%edge_index( vj) == 3 .OR. mesh%edge_index( vj) == 4)) THEN
        ! East
        mesh%edge_index_Ac( aci) = 3
      ELSEIF ((mesh%edge_index( vi)==4 .OR. mesh%edge_index( vi) == 5 .OR. mesh%edge_index( vi) == 6) .AND. &
              (mesh%edge_index( vj)==4 .OR. mesh%edge_index( vj) == 5 .OR. mesh%edge_index( vj) == 6)) THEN
        ! East
        mesh%edge_index_Ac( aci) = 5
      ELSEIF ((mesh%edge_index( vi)==6 .OR. mesh%edge_index( vi) == 7 .OR. mesh%edge_index( vi) == 8) .AND. &
              (mesh%edge_index( vj)==6 .OR. mesh%edge_index( vj) == 7 .OR. mesh%edge_index( vj) == 8)) THEN
        ! East
        mesh%edge_index_Ac( aci) = 7
      ELSE
        ! Not an edge Ac vertex
        mesh%edge_index_Ac( aci) = 0
      END IF
      
    END DO ! DO aci = mesh%ci1, mesh%ci2
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)
    
  END SUBROUTINE find_Ac_edge_indices

END MODULE mesh_ArakawaC_module
