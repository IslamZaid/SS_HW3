program Sphere_Model

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: numNodes = 16
  REAL(kind=8), PARAMETER :: r_ear = 6.371e+6 
  REAL(kind=8), PARAMETER :: pi = 3.1415926535897932388d0
  REAL(kind=8) :: X(numNodes, numNodes), Y(numNodes, numNodes), Z(numNodes, numNodes)
  REAL(kind=8) :: centerNodes(numNodes, numNodes, 4)
  REAL(kind=8) :: Area(numNodes, numNodes)
  INTEGER :: i, j
  REAL(kind=8) :: V1(3), V2(3), V3(3)
  REAL(kind=8) :: C1(3), C2(3) 
  REAL(kind=8) :: node1(3), node2(3), node3(3), node4(3)
  REAL(kind=8) :: AA
  real(kind=8) :: normalvector(numNodes, numNodes,3)
  real(kind=8) :: b_ear_sun, rad_ear_sun, rad_sat_ear,  area_sat, R_sat_sun(3), R_ear_sun(3), R_sat_ear(3)
  real(kind=8) :: dot_prod,  a, Ps, d
  real(kind=8) :: II(numNodes, numNodes), Q(numNodes, numNodes)

 
  ! Sphere Model Development
  X = 0.0d0
  Y = 0.0d0
  Z = 0.0d0


  DO i = 1, numNodes
    DO j = 1, numNodes
      X(i, j) = r_ear * SIN((i-1) * pi / (numNodes-1)) * COS((j-1) * 2 * pi / (numNodes-1))
      Y(i, j) = r_ear * SIN((i-1) * pi / (numNodes-1)) * SIN((j-1) * 2 * pi / (numNodes-1))
      Z(i, j) = r_ear * COS((i-1) * pi / (numNodes-1))
    END DO
  END DO


  centerNodes = 0.0d0
  Area = 0.0d0



  DO i = 1, numNodes-1
    DO j = 1, numNodes-1
      
      node1 = [X(i, j), Y(i, j), Z(i, j)]
      node2 = [X(i, j+1), Y(i, j+1), Z(i, j+1)]
      node3 = [X(i+1, j), Y(i+1, j), Z(i+1, j)]
      node4 = [X(i+1, j+1), Y(i+1, j+1), Z(i+1, j+1)]
      
      V1 = node2 - node1
      V2 = node3 - node1
      V3 = node4 - node1
      
      C1(1) = V1(2) * V2(3) - V1(3) * V2(2)
      C1(2) = V1(3) * V2(1) - V1(1) * V2(3)
      C1(3) = V1(1) * V2(2) - V1(2) * V2(1)
      
      
      C2(1) = V2(2) * V3(3) - V2(3) * V3(2)
      C2(2) = V2(3) * V3(1) - V2(1) * V3(3)
      C2(3) = V2(1) * V3(2) - V2(2) * V3(1)
      
      AA = NORM2(C1) + NORM2(C1)
      
      !print *, i,j,":", C1/Norm2(C1)
      
      Area(i, j) = AA
      
      ! Save center node information
      centerNodes(i, j, 1) = (X(i, j) + X(i, j+1) + X(i+1, j) + X(i+1, j+1)) / 4
      centerNodes(i, j, 2) = (Y(i, j) + Y(i, j+1) + Y(i+1, j) + Y(i+1, j+1)) / 4
      centerNodes(i, j, 3) = (Z(i, j) + Z(i, j+1) + Z(i+1, j) + Z(i+1, j+1)) / 4
      centerNodes(i, j, 4) = Area(i, j)
      
      ! Add label to center node
      
      normalvector(i, j, :) = C1/Norm2(C1)
      
    END DO
  END DO
  
  
  
  ! Radiation models
    ! position vectors definition
    b_ear_sun = pi / 2.d0
    rad_ear_sun = 1.496d11
    rad_sat_ear = 6.371d6 + 2.0d6
    area_sat = 1.d0
  
    R_ear_sun = rad_ear_sun * [cos(b_ear_sun), sin(b_ear_sun), 0.d0]
    R_sat_ear = rad_sat_ear * [cos(b_ear_sun), sin(b_ear_sun), 0.d0]
    R_sat_sun = R_ear_sun + R_sat_ear
    a = 0.39
    Ps = 1.00
  ! find the elements that see the sun and the radiant power emitted by sun per unit projected area
  DO i = 1 , numNodes
    DO j = 1 , numNodes
      V1 = normalvector(i,j,:)
      V2 = rad_ear_sun * [cos(b_ear_sun), sin(b_ear_sun), 0.d0] - [X(i, j), Y(i, j), Z(i, j)]
      dot_prod = V1(1)*V2(1) + V1(2)*V2(2) + V1(3)*V2(3) 
      d = norm2( rad_ear_sun * [cos(b_ear_sun), sin(b_ear_sun), 0.d0] - [X(i, j), Y(i, j), Z(i, j)] )
      if (dot_prod<0) then
        II(i,j) = Ps / ( 4*pi*d**2 )
      else 
        II(i,j) = 0
      end if
    END DO 
  END DO 
  
  ! Calculate the albedo 
  DO i = 1 , numNodes
    DO j = 1 , numNodes
      V1 = normalvector(i,j,:)
      V2 = rad_ear_sun * [cos(b_ear_sun), sin(b_ear_sun), 0.d0] - [X(i, j), Y(i, j), Z(i, j)]
      dot_prod = V1(1)*V2(1) + V1(2)*V2(2) + V1(3)*V2(3) 
      d = norm2( rad_ear_sun * [cos(b_ear_sun), sin(b_ear_sun), 0.d0] - [X(i, j), Y(i, j), Z(i, j)] )
      Q(i,j) = a * II(i,j) * dot_prod / d**2
	  
	  if Q(i,j) 
    END DO 
  END DO 
  
  
print *, normalvector(1,1,:)



end program Sphere_Model
