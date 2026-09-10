chrisharing@dhcp-168-105-243-106 hspist3 % 

26_10_02

tunneling at spring 2 wall: 

./00ALLINONE --mode=time --show-simulation --experiment=energy_transfer  --piston-right-protocol-mode=linear --max-right-piston-travel=10.0 --piston-right-duration=2.0 --velocity-right-piston-max=100.0 --velocity-right-piston-t0=0.0 --velocity-right-piston-gradient=10.0 --wall-hold-steps=10000 --steps=30000 --energy-measurement --spring-k=10.0 --spring-eq=20.0 --l0=20.0 --height=10.0 --output-dt=1.0 --left-empty --wall-thickness=0.05 --wall-thickness-vis=0.05 --wall-mass-factor=200.0 --kbt1 --eta=0.2 --particles=1000 --seed=2179731470 --num-walls=2 --wall-positions=20.0,30.0 --wall-mass-factors=200,200 --particles-boxes=0,500,500


--> change wall size

./00ALLINONE --mode=time --show-simulation --experiment=energy_transfer  --piston-right-protocol-mode=linear --max-right-piston-travel=10.0 --piston-right-duration=2.0 --velocity-right-piston-max=100.0 --velocity-right-piston-t0=0.0 --velocity-right-piston-gradient=10.0 --wall-hold-steps=10000 --steps=30000 --energy-measurement --spring-k=10.0 --spring-eq=20.0 --l0=20.0 --height=10.0 --output-dt=1.0 --left-empty --wall-thickness=0.4 --wall-thickness-vis=0.4  --wall-mass-factor=200.0 --kbt1 --eta=0.2 --particles=1000 --seed=2179731470 --num-walls=2 --wall-positions=20.0,30.0 --wall-mass-factors=200,200 --particles-boxes=0,500,500




1WALL:
./00ALLINONE --mode=time --show-simulation --experiment=energy_transfer  --piston-right-protocol-mode=linear --max-right-piston-travel=10.0 --piston-right-duration=2.0 --velocity-right-piston-max=100.0 --velocity-right-piston-t0=0.0 --velocity-right-piston-gradient=10.0 --wall-hold-steps=10000 --steps=30000 --energy-measurement --spring-k=10.0 --spring-eq=20.0 --l0=20.0 --height=10.0 --output-dt=1.0 --left-empty --wall-thickness=0.4 --wall-thickness-vis=0.4  --wall-mass-factor=200.0 --kbt1 --eta=0.2 --particles=1000 --seed=2179731470 --num-walls=1 --wall-positions=20.0 --wall-mass-factors=200 --particles-boxes=0,1000

./00ALLINONE --mode=time --show-simulation --experiment=energy_transfer  --piston-right-protocol-mode=linear --max-right-piston-travel=10.0 --piston-right-duration=2.0 --velocity-right-piston-max=100.0 --velocity-right-piston-t0=0.0 --velocity-right-piston-gradient=10.0 --wall-hold-steps=10000 --steps=30000 --energy-measurement --spring-k=20.0 --spring-eq=20.0 --l0=20.0 --height=10.0 --output-dt=1.0 --left-empty --wall-thickness=0.4 --wall-thickness-vis=0.4  --wall-mass-factor=200.0 --kbt1 --eta=0.5 --particles=1000 --seed=2179731470 --num-walls=1 --wall-positions=20.0 --wall-mass-factors=200 --particles-boxes=0,5000


few particles:

./00ALLINONE --mode=time --show-simulation --experiment=energy_transfer  --piston-right-protocol-mode=linear --max-right-piston-travel=10.0 --piston-right-duration=2.0 --velocity-right-piston-max=100.0 --velocity-right-piston-t0=0.0 --velocity-right-piston-gradient=10.0 --wall-hold-steps=10000 --steps=30000 --energy-measurement --spring-k=10.0 --spring-eq=20.0 --l0=20.0 --height=10.0 --output-dt=1.0 --left-empty --wall-thickness=0.4 --wall-thickness-vis=0.4  --wall-mass-factor=200.0 --kbt1 --seed=2179731470 --num-walls=2 --wall-positions=20.0,30.0 --wall-mass-factors=200,200 --particles-boxes=0,5,5 --particle-radius=1

./00ALLINONE --mode=time --show-simulation --experiment=energy_transfer  --piston-right-protocol-mode=linear --max-right-piston-travel=10.0 --piston-right-duration=2.0 --velocity-right-piston-max=100.0 --velocity-right-piston-t0=0.0 --velocity-right-piston-gradient=10.0 --wall-hold-steps=10000 --steps=30000 --energy-measurement --spring-k=10.0 --spring-eq=20.0 --l0=20.0 --height=10.0 --output-dt=1.0 --left-empty --wall-thickness=0.4 --wall-thickness-vis=0.4  --wall-mass-factor=200.0 --kbt1 --seed=2179731470 --num-walls=1 --wall-positions=20.0 --wall-mass-factors=200 --particles-boxes=0,10 --particle-radius=1


COMPLEX COUPLED - MULTIOPLE DIVIDER SPRINGS

2 WALL time
./00ALLINONE --mode=time\
  --num-walls=2 --wall-positions=20,30 --wall-mass-factors=200 \
  --energy-measurement --spring-k=10 --spring-eq=20 \
  --particles=100 --l0=20 --height=10 --kbt1 \
  --output-dt=1 --left-empty --particles-boxes=0,50,50 --wall-thickness=0.05 --wall-thickness-vis=0.05 --eta=0.196

--> around 90% energy in spring as in piston in max

1 Wall TIME

  chrisharing@dhcp-168-105-243-106 hspist3 % 
  
  ./00ALLINONE --mode=time\
  --num-walls=1 --wall-positions=20 --wall-mass-factors=200 \
  --energy-measurement --spring-k=10 --spring-eq=20 \
  --particles=100 --l0=20 --height=10 --kbt1 \
  --output-dt=1 --left-empty --particles-boxes=0,100 --wall-thickness=0.05 --wall-thickness-vis=0.05 --eta=0.196

  --> 80-90% piston in

chrisharing@dhcp-168-105-243-106 hspist3 % ./00ALLINONE --mode=time\
  --num-walls=3 --wall-positions=20,27,34 --wall-mass-factors=200 \
  --energy-measurement --spring-k=10 --spring-eq=20 \
  --particles=99 --l0=20 --height=10 --kbt1 \
  --output-dt=1 --left-empty --particles-boxes=0,33,33,33 --wall-thickness=0.1 --wall-thickness-vis=0.1 --eta=0.196 --seeding=random

  1 Wall edmd

  chrisharing@dhcp-168-105-243-106 hspist3 % ./00ALLINONE --mode=edmd\
  --num-walls=1 --wall-positions=20 --wall-mass-factors=200 \
  --energy-measurement --spring-k=10 --spring-eq=20 \
  --particles=100 --l0=20 --height=10 --kbt1 \
  --output-dt=1 --left-empty --particles-boxes=0,100 --wall-thickness=0.05 --wall-thickness-vis=0.05 --eta=0.196

  -->99.9 or more energy in spring than piston in





  /11/21/25

  chrisharing@dhcp-168-105-243-106 hspist3 % ./00ALLINONE --mode=time\
  --num-walls=3 --wall-positions=20,27,34 --wall-mass-factors=200 \
  --energy-measurement --spring-k=10 --spring-eq=20 \
  --particles=99 --l0=20 --height=10 --kbt1 \
  --output-dt=1 --left-empty --particles-boxes=0,33,33,33 --wall-thickness=1.0 --wall-thickness-vis=1.0 --eta=0.196 --seeding=random

  0.92


  chrisharing@dhcp-168-105-243-106 hspist3 % ./00ALLINONE --mode=time\
  --num-walls=2 --wall-positions=20,30 --wall-mass-factors=200 \
  --energy-measurement --spring-k=10 --spring-eq=20 \
  --particles=100 --l0=20 --height=10 --kbt1 \
  --output-dt=1 --left-empty --particles-boxes=0,50,50 --wall-thickness=1.0 --wall-thickness-vis=1.0 --eta=0.196 --seeding=random
  

  0.927




  SPEED OF SOUND - CS

   ./00ALLINONE --mode=edmd --particles=100 --l0=20 --height=10 \
  --wall-mass-factor=200 --temperature=100000 --kbt1 \
  --auto-release --steps=1500000 --single-test

  