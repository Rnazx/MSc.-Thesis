echo "Input the value of \Sigmatot"
read dat_sigmatot
echo "Input the value of \SigmaHI"

read dat_sigmaHI
echo "Input the value of q"
read dat_q
echo "Input the value of Omega"
read dat_omega
echo "Input the value of \sigmasfr"
read dat_sigmasfr
echo "Input the value of molfraction"
read molfrac
python arbi.py $dat_sigmatot $dat_sigmaHI $dat_q $dat_omega $dat_sigmasfr $molfrac
python arb_mag.py