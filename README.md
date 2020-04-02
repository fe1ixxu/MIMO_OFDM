# k-Modular Quadratic Programming Algorithm for PAPR in MIMO OFDM
we propose an alternative peak-toaveragepower ratio (PAPR) reduction framework for MIMO-OFDM
system based on the well-known unimodular quadratic
programming (UQP). In addition, we consider a more general
setting for PAPR reduction problem in MIMO-OFDM systems
and propose a novel power method-like algorithm to effectively
tackle the associated UQP. The proposed method can handle
arbitrary peak-to-average-power ratio (PAPR) constraints on
the transmit sequence, and more importantly, can be used to
generate constant modulus signals for such systems. The proposed
algorithm demonstrates an improvement in terms of convergence
rate compared with the state-of-the-art PAPR reduction method.

<div align=center>
  
<img src="https://github.com/fe1ixxu/MIMO_OFDM/blob/master/pictures/papr_10000.jpg" alt="PAPR" width="512px">
<img src="https://github.com/fe1ixxu/MIMO_OFDM/blob/master/pictures/papr_iteration.jpg" alt="one_block_iteartion" width="512px">
<img src="https://github.com/fe1ixxu/MIMO_OFDM/blob/master/pictures/BER.jpg" alt="BER" width="512px">
</div>

# Platform
MATLAB 2018b

# Simulation
Run the code and see the perfomance of k-MQP, run file

```
OFDM_Main.m
```

To see the BER, run file

```
BER.m
```
