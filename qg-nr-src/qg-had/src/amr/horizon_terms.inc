      integer soc(4)
        integer i40,sz
        parameter(i40=35,sz=2197)
        data soc/4,10,20,35/
        integer coeffsize
c			common/expol/soc
           parameter(coeffsize=20 000 000)
c         parameter(coeffsize=10700000)
c maximum of 2000 points at 4th order ....
c         parameter(coeffsize=6000)
c        parameter(coeffsize=700000)
