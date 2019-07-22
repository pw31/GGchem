      program JANAFFIT
      implicit none
      integer N,S
      parameter (N=200, S=10)
      integer specie,Nmax(S),Ndat,dN,i,j,grad,ed,Edzahl,mode,sumst
      integer Nit,it,NB
      real*8 dG(S,N),T(S,N),x(N),y(N),stoich(S),dGform(N),TTform(N)
      real*8 T1,T2,koeff(0:10),xx,yy,yfit,delta1,delta2
      real*8 bk,Rgas,bar,xlow,xhigh,Tdatmin,Tdatmax,Ggrad
      real*8 Afit(3),Bfit(5),e,ee1,dev
      real*8,external :: pvap,stock     ! function to fit
      logical gef,ok
      character(len=1) :: answer
      data bk/1.380662D-16/, Rgas/8.3144598D+0/, bar/1.D+6/
*
      write(*,*)
      write(*,*) 'Fit for which species?'
      write(*,*) '======================'
      write(*,*)
      write(*,*) ' 1 = CO '
      write(*,*) ' 2 = COS '
      write(*,*) ' 3 = Fe[s] '
      write(*,*) ' 4 = FeS[s] '
      write(*,*) ' 5 = Cl+ '
      write(*,*) ' 6 = Cl- '
      write(*,*) ' 7 = NaCl '
      write(*,*) ' 8 = NaCl[s] '
      write(*,*) ' 9 = Na[s] '
      write(*,*) '10 = CaCl '
      write(*,*) '11 = CaCl2 '
      write(*,*) '12 = KCl '
      write(*,*) '13 = LiH '
      write(*,*) '14 = LiO '
      write(*,*) '15 = LiOH '
      write(*,*) '16 = LiCl '
      write(*,*) '17 = HCl '
      write(*,*) '18 = MgS '
      write(*,*) '19 = MgS[s] '
      write(*,*) '20 = LiCl[s] '
      write(*,*) '21 = KCl[s]'
      write(*,*) '22 = Na2SiO3[s]'
      write(*,*) '23 = SiS2[s]'      
      write(*,*) '24 = CaS[s]'      
      write(*,*) '25 = FeS2[s]'      
      write(*,*) '26 = H2O[l]'      
      write(*,*) '27 = H2SO4'      
      write(*,*) '28 = H2SO4[l]'      
      write(*,*) '29 = C4N2'      
      write(*,*) '30 = C[s]'      
      write(*,*) '31 = CH4'      
      write(*,*) '32 = Na2S[s]'      
      write(*,*) '33 = MgCl'      
      write(*,*) '34 = MgCl2'
      write(*,*) '35 = Mg2Cl4'
      write(*,*) '36 = TiCl2'
      write(*,*) '37 = TiCl4'
      write(*,*) '38 = FeCl'
      write(*,*) '39 = FeCl2'
      write(*,*) '40 = AlCl'
      write(*,*) '41 = AlCl2'
      write(*,*) '42 = AlCl3'
      write(*,*) '43 = Al2Cl6'
      write(*,*) '44 = SiCl'
      write(*,*) '45 = SiCl2'
      write(*,*) '46 = SiCl3'
      write(*,*) '47 = SiCl4'
      write(*,*) '48 = N+'
      write(*,*) '49 = AlCl3[s]'
      write(*,*) '50 = HCO'
      write(*,*) '51 = NH2'
      write(*,*) '52 = NH'
      write(*,*) '53 = CO2'
      write(*,*) '54 = HCN'
      write(*,*) '55 = N2'
      write(*,*) '56 = NH3'
      write(*,*) '57 = CH2'
      write(*,*) '58 = CH3'
      write(*,*) '59 = H2'
      write(*,*) '60 = H-'
      write(*,*) '61 = H2-'
      write(*,*) '62 = H2+'
      write(*,*) '63 = He+'
      write(*,*) '64 = SiN'
      write(*,*) '65 = CS'
      write(*,*) '66 = SiS'
      write(*,*) '67 = HS does not work!'
      write(*,*) '68 = H2S'
      write(*,*) '69 = SiO'
      write(*,*) '70 = SiC2'
      write(*,*) '71 = Na2O2H2'
      write(*,*) '72 = NaOH'
      write(*,*) '73 = NaCN'
      write(*,*) '74 = Na2'
      write(*,*) '75 = S2'
      write(*,*) '76 = C2H4'
      write(*,*) '77 = SN'
      write(*,*) '78 = CN'
      write(*,*) '79 = K2O2H2'
      write(*,*) '80 = KOH'
      write(*,*) '81 = AlOH'
      write(*,*) '82 = HAlO'
      write(*,*) '83 = Al2O'
      write(*,*) '84 = AlO2H'
      write(*,*) '85 = Al2O2'
      write(*,*) '86 = AlS'
      write(*,*) '87 = AlO'
      write(*,*) '88 = MgOH'
      write(*,*) '89 = MgO2H2'
      write(*,*) '90 = MgO'
      write(*,*) '91 = FeS'
      write(*,*) '92 = AlH'
      write(*,*) '93 = H2O'
      write(*,*) '94 = SiH4'
      write(*,*) '95 = TiO2'
      write(*,*) '96 = TiO'
      write(*,*) '97 = Ca(OH)2'
      write(*,*) '98 = K+'
      write(*,*) '99 = Na+'
      write(*,*) '100 = Ca+'
      write(*,*) '101 = AlO2-'
      write(*,*) '102 = CN-'
      write(*,*) '103 = Li+'
      write(*,*) '104 = AlOH-'
      write(*,*) '105 = CaO'
      write(*,*) '106 = CKN'
      write(*,*) '107 = CaS'
      write(*,*) '108 = AlOH+'
      write(*,*) '109 = S8'
      write(*,*) '110 = H3O+'
      write(*,*) '111 = Al2O+'
      write(*,*) '112 = MgOH++'
      write(*,*) '113 = N2H4'
      write(*,*) '114 = Mg+'      
      write(*,*) '115 = Al+'      
      write(*,*) '116 = Al2O2+'
      write(*,*) '117 = AlO+'
      write(*,*) '118 = HCO+'
      write(*,*) '119 = NO+'
      write(*,*) '120 = H+'
      write(*,*) '121 = FeO2H2'
      write(*,*) '122 = SiO2'
      write(*,*) '123 = OH-'
      write(*,*) '124 = CO2-'
      write(*,*) '125 = AlO-'
      write(*,*) '126 = NaO-'
      write(*,*) '127 = MgO[s]'
      write(*,*) '128 = C2H4O'
      write(*,*) '129 = SiH'
      write(*,*) '130 = Li2O2H2'
      write(*,*) '131 = Li2Cl2'
      write(*,*) '132 = Li3Cl3'
      write(*,*) '133 = Fe(CO)5'
      write(*,*) '134 = Si(CH3)4'
      write(*,*) '135 = (KCl)2'
      write(*,*) '136 = (NaCl)2'
      write(*,*) '137 = TiOCl2'
      write(*,*) '138 = TiOCl'
      write(*,*) '139 = Na2SO4'
      write(*,*) '140 = K2SO4'
      write(*,*) '141 = S7'
      write(*,*) '142 = AlOCl'
      write(*,*) '143 = KO'
      write(*,*) '144 = LiN'
      write(*,*) '145 = SiCH3Cl3'
      write(*,*) '146 = Fe[l]'
      write(*,*) '147 = FeS[l]'
      write(*,*) '148 = NaCl[l]'
      write(*,*) '149 = Na[l]'
      write(*,*) '150 = CaCl2[l]'
      write(*,*) '151 = KCl[l]'
      write(*,*) '152 = LiH[l]'
      write(*,*) '153 = LiOH[l]'
      write(*,*) '154 = LiCl[l]'
      write(*,*) '155 = Na2SiO3[l]'
      write(*,*) '156 = SiS2[l]'
      write(*,*) '157 = Na2S[l]'
      write(*,*) '158 = MgCl2[l]'
      write(*,*) '159 = TiCl4[l]'
      write(*,*) '160 = FeCl2[l]'
      write(*,*) '161 = AlCl3[l]'
      write(*,*) '162 = NaOH[l]'
      write(*,*) '163 = NaCN[l]'
      write(*,*) '164 = KOH[l]'
      write(*,*) '165 = MgO[l]'
      write(*,*) '166 = TiO2[l]'
      write(*,*) '167 = TiO[l]'
      write(*,*) '168 = CaO[l]'
      write(*,*) '169 = CKN[l]'
      write(*,*) '170 = N2H4[l]'
      write(*,*) '171 = SiO2[l]'
      write(*,*) '172 = TiO2[s]'
      write(*,*) '173 = Fe(CO)5[l]'
      write(*,*) '174 = Na2SO4[l]'
      write(*,*) '175 = K2SO4[l]'
      write(*,*) '176 = Al2O3[s]'
      write(*,*) '177 = Al2O3[l]'
      write(*,*) '178 = MgAl2O4[s]'
      write(*,*) '179 = MgAl2O4[l]'
      write(*,*) '180 = Ti4O7[s]'
      write(*,*) '181 = Ti4O7[l]'
      write(*,*) '182 = Mg2SiO4[l]'
      write(*,*) '183 = MgSiO3[l]'
      write(*,*) '184 = CaO[s]'
      write(*,*) '185 = TiO2[s] pvap'
      write(*,*) '186 = TiO2[l] pvap'
      write(*,*) '187 = H2O[l] pvap'
      write(*,*) '188 = SiH3Cl'
      write(*,*) '189 = FeCl3'
      write(*,*) '190 = MgSiO3[s]'
      write(*,*) '191 = SiO2[s]'
      write(*,*) '192 = Mg2SiO4[s]'
      write(*,*) '193 = AlO2'
      write(*,*) '194 = C5'
      write(*,*) '195 = S2O'
      write(*,*) '196 = FeO[l]'
      write(*,*) '197 = FeO[s]'
      write(*,*) '198 = SiO2[s] pvap'
      write(*,*) '199 = SiO2[l] pvap'
      write(*,*) '200 = FeS[s] pvap'
      write(*,*) '201 = FeS[l] pvap'
      write(*,*) '202 = NaCl[s] pvap'
      write(*,*) '203 = NaCl[l] pvap'
      write(*,*) '204 = KCl[s] pvap'
      write(*,*) '205 = KCl[l] pvap'
      write(*,*) '206 = LiH[s] pvap'
      write(*,*) '207 = LiH[l] pvap'
      write(*,*) '208 = LiH[s]'
      write(*,*) '209 = CaCl2[s] pvap'
      write(*,*) '210 = CaCl2[l] pvap'
      write(*,*) '211 = CaCl2[s]'
      write(*,*) '212 = FeO'
      write(*,*) '213 = FeO[s] pvap'
      write(*,*) '214 = FeO[l] pvap'
      write(*,*) '215 = MgTi2O5[s]'
      write(*,*) '216 = MgTi2O5[l]'
      write(*,*) '217 = MgO[s] pvap'
      write(*,*) '218 = MgO[l] pvap'
      write(*,*) '219 = AlCl3[s] pvap'
      write(*,*) '220 = AlCl3[l] pvap'
      write(*,*) '221 = LiCl[s] pvap'
      write(*,*) '222 = LiCl[l] pvap'
      write(*,*) '223 = MgTiO3[s]'
      write(*,*) '224 = MgTiO3[l]'
      write(*,*) '225 = CaO[s] pvap'
      write(*,*) '226 = CaO[l] pvap'
      write(*,*) '227 = S[s]'
      write(*,*) '228 = S[l]'
      write(*,*) '229 = K2SiO3[s]'
      write(*,*) '230 = K2SiO3[l]'
      write(*,*) '231 = TiC[s]'
      write(*,*) '232 = TiC[l]'
      write(*,*) '233 = Ti[s]'
      write(*,*) '234 = Ti[l]'
      write(*,*) '235 = TiO[s]'
      write(*,*) '236 = TiO[s] pvap'
      write(*,*) '237 = TiO[l] pvap'
      write(*,*) '238 = LiOH[s]'
      write(*,*) '239 = LiOH[s] pvap'
      write(*,*) '240 = LiOH[l] pvap'
      write(*,*) '241 = TiCl3'
      write(*,*) '242 = atomic Al'
      write(*,*) '243 = atomic O'
      write(*,*) '244 = atomic Mg'
      write(*,*) '245 = atomic Ca'
      write(*,*) '246 = atomic Si'
      write(*,*) '247 = atomic Na'
      write(*,*) '248 = atomic K'
      write(*,*) '249 = atomic Ti'
      write(*,*) '250 = atomic Fe'
      write(*,*) '251 = atomic S'
      write(*,*) '252 = atomic C'
      write(*,*) '253 = atomic Cl'
      write(*,*) '254 = atomic Li'
      write(*,*) '255 = atomic H'
      write(*,*) '256 = atomic N'
      write(*,*) '272 = atomic Mn'
      write(*,*) '273 = atomic Cr'
      write(*,*) '274 = atomic Zr'
      write(*,*) '275 = atomic Cu'
      write(*,*) '276 = atomic F'
      write(*,*) '277 = atomic P'
      write(*,*) '278 = atomic Ni'
      write(*,*) '257 = W[s] pvap'
      write(*,*) '258 = W[l] pvap'
      write(*,*) '259 = WO3[s] pvap'
      write(*,*) '260 = WO3[l] pvap'
      write(*,*) '261 = ZrO2[s] pvap'
      write(*,*) '262 = ZrO2[l] pvap'
      write(*,*) '263 = VO[s] pvap'
      write(*,*) '264 = VO[l] pvap'
      write(*,*) '265 = Zr[s] pvap'
      write(*,*) '266 = Zr[l] pvap'
      write(*,*) '267 = ZrO2[s]'
      write(*,*) '268 = ZrO2[l]'
      write(*,*) '269 = SiC[s]'
      write(*,*) '270 = Mn[s]'
      write(*,*) '271 = Mn[l]'
      write(*,*) '279 = PH3'
      write(*,*) '280 = PH2'
      write(*,*) '281 = PH'
      write(*,*) '282 = PN'
      write(*,*) '283 = SO2'
      write(*,*) '284 = MgClF'
      write(*,*) '285 = WO'
      write(*,*) '286 = WO2'
      write(*,*) '287 = WO3'
      write(*,*) '288 = WCl2'
      write(*,*) '289 = ZrSiO4_cr'
      write(*,*) '290 = H2WO4'
      write(*,*) '291 = WCl'
      write(*,*) '292 = WF'
      write(*,*) '293 = WO2Cl2'
      write(*,*) '294 = W2O6'
      write(*,*) '295 = W3O8'
      write(*,*) '296 = W3O9'
      write(*,*) '297 = W4O12'
      write(*,*) '298 = Ni_cr'
      write(*,*) '299 = Ni_l'
      write(*,*) '300 = Cr_cr'
      write(*,*) '301 = Cr_l'
      write(*,*) '302 = CrN_cr pvap'
      write(*,*) '303 = V2O3_cr'
      write(*,*) '304 = V2O4_cr'
      write(*,*) '305 = V2O5_cr'
      write(*,*) '306 = Ni3S2_cr'
      write(*,*) '307 = Ni3S2_l'
      write(*,*) '308 = NH4Cl_cr'
      write(*,*) '309 = P_cr'
      write(*,*) '310 = P_l'
      write(*,*) '311 = P2O10_cr'
      write(*,*) '312 = P4S3_cr'
      write(*,*) '313 = P4S3_l'
      write(*,*) '314 = VO'
      write(*,*) '315 = VO2'
      write(*,*) '316 = SO3'
      write(*,*) '317 = Zn_cr'
      write(*,*) '318 = Zn_l'
      write(*,*) '319 = ZnSO4_cr'
      write(*,*) '320 = H3PO4_cr'
      write(*,*) '321 = Mg3P2O8_cr'
      write(*,*) '322 = P3N5_cr'
      write(*,*) '323 = P4O6'
      write(*,*) '324 = P4O10'
      write(*,*) '325 = AlF3_cr'
      write(*,*) '326 = CaF2 pvap'
      write(*,*) '327 = KF pvap'
      write(*,*) '328 = NaF pvap'
      write(*,*) '329 = H3PO4_l'
      write(*,*) '330 = FeF2 pvap'
      write(*,*) '331 = MgF2 pvap'      
      write(*,*) '332 = AlF6Na3_cr'      
      write(*,*) '333 = Li2SiO3_cr'      
      write(*,*) '334 = Li2SiO3_l'      
      write(*,*) '335 = Li2Si2O5_cr'      
      write(*,*) '336 = Li2Si2O5_l'      
      write(*,*) '337 = Li2TiO3_cr'      
      write(*,*) '338 = Li2TiO3_l'      
      write(*,*) '339 = Co_cr'
      write(*,*) '340 = Co_l'
      write(*,*) '341 = CoO_cr'
      read(*,*) specie
*
      if (specie.eq.1) then
        call READ_DATEI('CO.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.2) then
        call READ_DATEI('COS.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.3) then
        call READ_DATEI('Fe_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.4) then
        call READ_DATEI('FeS_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.5) then
        call READ_DATEI('Cl+.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Cl.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt',dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = -1.D0
      elseif (specie.eq.6) then
        call READ_DATEI('Cl-.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Cl.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt',dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.7) then
        call READ_DATEI('NaCl.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Cl.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Na.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.8) then
        call READ_DATEI('NaCl_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Cl.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Na.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.9) then
        call READ_DATEI('Na_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.10) then
        call READ_DATEI('CaCl.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.11) then
        call READ_DATEI('CaCl2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.12) then
        call READ_DATEI('KCl.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.13) then
        call READ_DATEI('LiH.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.14) then
        call READ_DATEI('LiO.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.15) then
        call READ_DATEI('LiOH.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.16) then
        call READ_DATEI('LiCl.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.17) then
        call READ_DATEI('HCl.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.18) then
        call READ_DATEI('MgS.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.19) then
        call READ_DATEI('MgS_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.20) then
        call READ_DATEI('LiCl_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.21) then
        call READ_DATEI('KCl_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.22) then
        call READ_DATEI('Na2SiO3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'     ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'      ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.23) then
        call READ_DATEI('SiS2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'      ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.24) then
        call READ_DATEI('CaS_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.25) then
        call READ_DATEI('FeS2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'      ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.26) then
        call READ_DATEI('H2O_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 1.D0
      elseif (specie.eq.27) then
        call READ_DATEI('H2SO4.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.28) then
        call READ_DATEI('H2SO4_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.29) then
        call READ_DATEI('C4N2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 4.D0
        stoich(3) = 2.D0
      elseif (specie.eq.30) then
        call READ_DATEI('C_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.31) then
        call READ_DATEI('CH4.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 4.D0
      elseif (specie.eq.32) then
        call READ_DATEI('Na2S_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'      ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 1.D0
      elseif (specie.eq.33) then
        call READ_DATEI('MgCl.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.34) then
        call READ_DATEI('MgCl2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.35) then
        call READ_DATEI('Mg2Cl4.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 4.D0
      elseif (specie.eq.36) then
        call READ_DATEI('TiCl2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.37) then
        call READ_DATEI('TiCl4.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 4.D0
      elseif (specie.eq.38) then
        call READ_DATEI('FeCl.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.39) then
        call READ_DATEI('FeCl2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.40) then
        call READ_DATEI('AlCl.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.41) then
        call READ_DATEI('AlCl2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.42) then
        call READ_DATEI('AlCl3.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.43) then
        call READ_DATEI('Al2Cl6.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 6.D0
      elseif (specie.eq.44) then
        call READ_DATEI('SiCl.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.45) then
        call READ_DATEI('SiCl2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.46) then
        call READ_DATEI('SiCl3.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.47) then
        call READ_DATEI('SiCl4.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 4.D0
      elseif (specie.eq.48) then
        call READ_DATEI('N+.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('N.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt',dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = -1.D0
      elseif (specie.eq.49) then
        call READ_DATEI('AlCl3_cr.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.50) then
        call READ_DATEI('HCO.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.51) then
        call READ_DATEI('NH2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.52) then
        call READ_DATEI('NH.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.53) then
        call READ_DATEI('CO2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.54) then
        call READ_DATEI('HCN.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.55) then
        call READ_DATEI('N2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 2.D0
      elseif (specie.eq.56) then
        call READ_DATEI('NH3.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.57) then
        call READ_DATEI('CH2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.58) then
        call READ_DATEI('CH3.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.59) then
        call READ_DATEI('H2.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 2.D0
      elseif (specie.eq.60) then
        call READ_DATEI('H-.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.61) then
        call READ_DATEI('H2-.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 1.D0
      elseif (specie.eq.62) then
        call READ_DATEI('H2+.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = -1.D0
      elseif (specie.eq.63) then
        call READ_DATEI('He+.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('He.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = -1.D0
      elseif (specie.eq.64) then
        call READ_DATEI('SiN.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.65) then
        call READ_DATEI('CS.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.66) then
        call READ_DATEI('SiS.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.67) then
        call READ_DATEI('HScorr.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.68) then
        call READ_DATEI('H2S.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 1.D0
      elseif (specie.eq.69) then
        call READ_DATEI('SiO.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.70) then
        call READ_DATEI('SiC2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.71) then
        call READ_DATEI('Na2O2H2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 2.D0
        stoich(4) = 2.D0
      elseif (specie.eq.72) then
        call READ_DATEI('NaOH.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.73) then
        call READ_DATEI('NaCN.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.74) then
        call READ_DATEI('Na2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'  ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 2.D0
      elseif (specie.eq.75) then
        call READ_DATEI('S2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 2.D0
      elseif (specie.eq.76) then
        call READ_DATEI('C2H4.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 4.D0
      elseif (specie.eq.77) then
        call READ_DATEI('SN.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('N.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.78) then
        call READ_DATEI('CN.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('N.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.79) then
        call READ_DATEI('K2O2H2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 2.D0
        stoich(4) = 2.D0
      elseif (specie.eq.80) then
        call READ_DATEI('KOH.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.81) then
        call READ_DATEI('AlOH.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.82) then
        call READ_DATEI('HAlO.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.83) then
        call READ_DATEI('Al2O.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 1.D0
      elseif (specie.eq.84) then
        call READ_DATEI('AlO2H.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 1.D0
      elseif (specie.eq.85) then
        call READ_DATEI('Al2O2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 2.D0
      elseif (specie.eq.86) then
        call READ_DATEI('AlS.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.87) then
        call READ_DATEI('AlO.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.88) then
        call READ_DATEI('MgOH.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.89) then
        call READ_DATEI('MgO2H2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 2.D0
      elseif (specie.eq.90) then
        call READ_DATEI('MgO.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.91) then
        call READ_DATEI('FeS.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.92) then
        call READ_DATEI('AlH.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.93) then
        call READ_DATEI('H2O.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 1.D0
      elseif (specie.eq.94) then
        call READ_DATEI('SiH4.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 4.D0
      elseif (specie.eq.95) then
        call READ_DATEI('TiO2.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.96) then
        call READ_DATEI('TiO.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.97) then
        call READ_DATEI('CaO2H2.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 2.D0
      elseif (specie.eq.98) then
        call READ_DATEI('K+.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = -1.D0
      elseif (specie.eq.99) then
        call READ_DATEI('Na+.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = -1.D0
      elseif (specie.eq.100) then
        call READ_DATEI('Ca+.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = -1.D0
      elseif (specie.eq.101) then
        call READ_DATEI('AlO2-.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 1.D0
      elseif (specie.eq.102) then
        call READ_DATEI('CN-.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('N.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.103) then
        call READ_DATEI('Li+.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = -1.D0
      elseif (specie.eq.104) then
        call READ_DATEI('AlOH-.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,4) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,5) 
        Edzahl = 4
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
        stoich(5) = 1.D0
      elseif (specie.eq.105) then
        call READ_DATEI('CaO.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.106) then
        call READ_DATEI('CKN.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('K.txt'  ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('N.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.107) then
        call READ_DATEI('CaS.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.108) then
        call READ_DATEI('AlOH+.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,4) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,5) 
        Edzahl = 4
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
        stoich(5) = -1.D0
      elseif (specie.eq.109) then
        call READ_DATEI('S8.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('S.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 8.D0
      elseif (specie.eq.110) then
        call READ_DATEI('H3O+.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 3.D0
        stoich(3) = 1.D0
        stoich(4) = -1.D0
      elseif (specie.eq.111) then
        call READ_DATEI('Al2O+.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = -1.D0
      elseif (specie.eq.112) then
        call READ_DATEI('MgOH+.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3)  
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,4) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,5) 
        Edzahl = 4
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
        stoich(5) = -1.D0
      elseif (specie.eq.113) then
        call READ_DATEI('N2H4.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('N.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,3)  
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 4.D0
      elseif (specie.eq.114) then
        call READ_DATEI('Mg+.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,3)  
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = -1.D0
      elseif (specie.eq.115) then
        call READ_DATEI('Al+.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,3)  
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = -1.D0
      elseif (specie.eq.116) then
        call READ_DATEI('Al2O2+.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'     ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 2.D0
        stoich(4) = -1.D0
      elseif (specie.eq.117) then
        call READ_DATEI('AlO+.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'     ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = -1.D0
      elseif (specie.eq.118) then
        call READ_DATEI('HCO+.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'     ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'     ,dG,T,Nmax,N,S,4) 
        call READ_DATEI('el-.txt'   ,dG,T,Nmax,N,S,5) 
        Edzahl = 4
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
        stoich(5) = -1.D0
      elseif (specie.eq.119) then
        call READ_DATEI('NO+.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('N.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = -1.D0
      elseif (specie.eq.120) then
        call READ_DATEI('H+.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = -1.D0
      elseif (specie.eq.121) then
        call READ_DATEI('FeO2H2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'     ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'     ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 2.D0
      elseif (specie.eq.122) then
        call READ_DATEI('SiO2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.123) then
        call READ_DATEI('OH-.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.124) then
        call READ_DATEI('CO2-.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt' ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 1.D0
      elseif (specie.eq.125) then
        call READ_DATEI('AlO-.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.126) then
        call READ_DATEI('NaO-.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('el-.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.127) then
        call READ_DATEI('MgO_cr.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.128) then
        call READ_DATEI('C2H4O.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 4.D0
        stoich(4) = 1.D0
      elseif (specie.eq.129) then
        call READ_DATEI('SiH.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.130) then
        call READ_DATEI('Li2O2H2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'      ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'      ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 2.D0
        stoich(4) = 2.D0
      elseif (specie.eq.131) then
        call READ_DATEI('Li2Cl2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 2.D0
      elseif (specie.eq.132) then
        call READ_DATEI('Li3Cl3.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 3.D0
        stoich(3) = 3.D0
      elseif (specie.eq.133) then
        call READ_DATEI('Fe(CO)5.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'      ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'      ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 5.D0
        stoich(4) = 5.D0
      elseif (specie.eq.134) then
        call READ_DATEI('Si(CH3)4.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'      ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'       ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'       ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 4.D0
        stoich(4) = 12.D0
      elseif (specie.eq.135) then
        call READ_DATEI('K2Cl2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 2.D0
      elseif (specie.eq.136) then
        call READ_DATEI('Na2Cl2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Cl.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Na.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 2.D0
      elseif (specie.eq.137) then
        call READ_DATEI('TiOCl2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('Cl.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 2.D0
      elseif (specie.eq.138) then
        call READ_DATEI('TiOCl.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('Cl.txt'  ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.139) then
        call READ_DATEI('Na2SO4.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.140) then
        call READ_DATEI('K2SO4.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.141) then
        call READ_DATEI('S7.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('S.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 7.D0
      elseif (specie.eq.142) then
        call READ_DATEI('AlOCl.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.143) then
        call READ_DATEI('KO.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.144) then
        call READ_DATEI('LiN.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('N.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.145) then
        call READ_DATEI('SiCH3Cl3.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,4) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,5) 
        Edzahl = 4
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
        stoich(5) = 3.D0
      elseif (specie.eq.146) then
        call READ_DATEI('Fe_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.147) then
        call READ_DATEI('FeS_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.148) then
        call READ_DATEI('NaCl_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.149) then
        call READ_DATEI('Na_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.150) then
        call READ_DATEI('CaCl2_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.151) then
        call READ_DATEI('KCl_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.152) then
        call READ_DATEI('LiH_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.153) then
        call READ_DATEI('LiOH_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.154) then
        call READ_DATEI('LiCl_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.155) then
        call READ_DATEI('Na2SiO3_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.156) then
        call READ_DATEI('SiS2_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.157) then
        call READ_DATEI('Na2S_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 1.D0
      elseif (specie.eq.158) then
        call READ_DATEI('MgCl2_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.159) then
        call READ_DATEI('TiCl4_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 4.D0
      elseif (specie.eq.160) then
        call READ_DATEI('FeCl2_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.161) then
        call READ_DATEI('AlCl3_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.162) then
        call READ_DATEI('NaOH_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.163) then
        call READ_DATEI('NaCN_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.164) then
        call READ_DATEI('KOH_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.165) then
        call READ_DATEI('MgO_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.166) then
        call READ_DATEI('TiO2_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.167) then
        call READ_DATEI('TiO_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.168) then
        call READ_DATEI('CaO_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.169) then
        call READ_DATEI('CKN_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('K.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.170) then
        call READ_DATEI('N2H4_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 4.D0
      elseif (specie.eq.171) then
        call READ_DATEI('SiO2_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.172) then
        call READ_DATEI('TiO2_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.173) then
        call READ_DATEI('Fe(CO)5_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 5.D0
        stoich(4) = 5.D0
      elseif (specie.eq.174) then
        call READ_DATEI('Na2SO4_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Na.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.175) then
        call READ_DATEI('K2SO4_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.176) then
        call READ_DATEI('Al2O3_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 3.D0
      elseif (specie.eq.177) then
        call READ_DATEI('Al2O3_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 3.D0
      elseif (specie.eq.178) then
        call READ_DATEI('MgAl2O4_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Al.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 4.D0
      elseif (specie.eq.179) then
        call READ_DATEI('MgAl2O4_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Al.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 4.D0
      elseif (specie.eq.180) then
        call READ_DATEI('Ti4O7_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 4.D0
        stoich(3) = 7.D0
      elseif (specie.eq.181) then
        call READ_DATEI('Ti4O7_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 4.D0
        stoich(3) = 7.D0
      elseif (specie.eq.182) then
        call READ_DATEI('Mg2SiO4_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.183) then
        call READ_DATEI('MgSiO3_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.184) then
        call READ_DATEI('CaO_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.185) then
        call READ_DATEI('TiO2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('TiO2.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.186) then
        call READ_DATEI('TiO2_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('TiO2.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.187) then
        call READ_DATEI('H2O_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('H2O.txt'  ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.188) then
        call READ_DATEI('SiH3Cl.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'     ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('Cl.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 3.D0
        stoich(4) = 1.D0
      elseif (specie.eq.189) then
        call READ_DATEI('FeCl3.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.190) then
        call READ_DATEI('MgSiO3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.191) then
        call READ_DATEI('SiO2_cr.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.192) then
        call READ_DATEI('Mg2SiO4_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.193) then
        call READ_DATEI('AlO2.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.194) then
        call READ_DATEI('C5.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('C.txt'  ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 5.D0
      elseif (specie.eq.195) then
        call READ_DATEI('S2Ocorr.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('S.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 1.D0
      elseif (specie.eq.196) then
        call READ_DATEI('FeO_l.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.197) then
        call READ_DATEI('FeO_cr.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'  ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'  ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.198) then
        call READ_DATEI('SiO2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('SiO2.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.199) then
        call READ_DATEI('SiO2_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('SiO2.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.200) then
        call READ_DATEI('FeS_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('FeS.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.201) then
        call READ_DATEI('FeS_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('FeS.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.202) then
        call READ_DATEI('NaCl_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('NaCl.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.203) then
        call READ_DATEI('NaCl_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('NaCl.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.204) then
        call READ_DATEI('KCl_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('KCl.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.205) then
        call READ_DATEI('KCl_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('KCl.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.206) then
        call READ_DATEI('LiH_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('LiH.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.207) then
        call READ_DATEI('LiH_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('LiH.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.208) then
        call READ_DATEI('LiH_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.209) then
        call READ_DATEI('CaCl2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('CaCl2.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.210) then
        call READ_DATEI('CaCl2_l.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('CaCl2.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.211) then
        call READ_DATEI('CaCl2_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ca.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.212) then
        call READ_DATEI('FeO.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Fe.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.213) then
        call READ_DATEI('FeO_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('FeO.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.214) then
        call READ_DATEI('FeO_l.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('FeO.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.215) then
        call READ_DATEI('MgTi2O5_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Ti.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 5.D0
      elseif (specie.eq.216) then
        call READ_DATEI('MgTi2O5_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Ti.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 5.D0
      elseif (specie.eq.217) then
        call READ_DATEI('MgO_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('MgO.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.218) then
        call READ_DATEI('MgO_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('MgO.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.219) then
        call READ_DATEI('AlCl3_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('AlCl3.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.220) then
        call READ_DATEI('AlCl3_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('AlCl3.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.221) then
        call READ_DATEI('LiCl_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('LiCl.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.222) then
        call READ_DATEI('LiCl_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('LiCl.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.223) then
        call READ_DATEI('MgTiO3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Ti.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.224) then
        call READ_DATEI('MgTiO3_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Ti.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.225) then
        call READ_DATEI('CaO_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('CaO.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.226) then
        call READ_DATEI('CaO_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('CaO.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.227) then
        call READ_DATEI('S_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('S.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.228) then
        call READ_DATEI('S_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('S.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.229) then
        call READ_DATEI('K2SiO3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.230) then
        call READ_DATEI('K2SiO3_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('K.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.231) then
        call READ_DATEI('TiC_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.232) then
        call READ_DATEI('TiC_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.233) then
        call READ_DATEI('Ti_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.234) then
        call READ_DATEI('Ti_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.235) then
        call READ_DATEI('TiO_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.236) then
        call READ_DATEI('TiO_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('TiO.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.237) then
        call READ_DATEI('TiO_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('TiO.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.238) then
        call READ_DATEI('LiOH_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3)
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.239) then
        call READ_DATEI('LiOH_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('LiOH.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.240) then
        call READ_DATEI('LiOH_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('LiOH.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.241) then
        call READ_DATEI('TiCl3.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.242) then
        call READ_DATEI('Alatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.243) then
        call READ_DATEI('Oatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.244) then
        call READ_DATEI('Mgatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.245) then
        call READ_DATEI('Caatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.246) then
        call READ_DATEI('Siatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.247) then
        call READ_DATEI('Naatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.248) then
        call READ_DATEI('Katom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.249) then
        call READ_DATEI('Tiatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.250) then
        call READ_DATEI('Featom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.251) then
        call READ_DATEI('Satom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.252) then
        call READ_DATEI('Catom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.253) then
        call READ_DATEI('Clatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.254) then
        call READ_DATEI('Liatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.255) then
        call READ_DATEI('Hatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.256) then
        call READ_DATEI('Natom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.272) then
        call READ_DATEI('Mnatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.273) then
        call READ_DATEI('Cratom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.274) then
        call READ_DATEI('Zratom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.275) then
        call READ_DATEI('Cuatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.276) then
        call READ_DATEI('Fatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.277) then
        call READ_DATEI('Patom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.278) then
        call READ_DATEI('Niatom.txt',dG,T,Nmax,N,S,1) 
        Edzahl = 0
      elseif (specie.eq.258) then
        call READ_DATEI('W_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt',dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.258) then
        call READ_DATEI('W_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt',dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.259) then
        call READ_DATEI('WO3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('WO3.txt',dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.260) then
        call READ_DATEI('WO3_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('WO3.txt',dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.261) then
        call READ_DATEI('ZrO2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('ZrO2.txt',dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.262) then
        call READ_DATEI('ZrO2_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('ZrO2.txt',dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.263) then
        call READ_DATEI('VO_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('VO.txt',dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.264) then
        call READ_DATEI('VO_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('VO.txt',dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.265) then
        call READ_DATEI('Zr_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Zr.txt',dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.266) then
        call READ_DATEI('Zr_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Zr.txt',dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.267) then
        call READ_DATEI('ZrO2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Zr.txt',dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt',dG,T,Nmax,N,S,3)
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.268) then
        call READ_DATEI('ZrO2_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Zr.txt',dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt',dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.269) then
        call READ_DATEI('SiC_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('C.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.270) then
        call READ_DATEI('Mn_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mn.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.271) then
        call READ_DATEI('Mn_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mn.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.279) then
        call READ_DATEI('PH3.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.280) then
        call READ_DATEI('PH2corr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.281) then
        call READ_DATEI('PH.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.282) then
        call READ_DATEI('PN.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('N.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.283) then
        call READ_DATEI('SO2.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('S.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.284) then
        call READ_DATEI('MgClF.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('F.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 1.D0
      elseif (specie.eq.285) then
        call READ_DATEI('WO.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.286) then
        call READ_DATEI('WO2.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.287) then
        call READ_DATEI('WO3.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.288) then
        call READ_DATEI('WCl2.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.289) then
        call READ_DATEI('ZrSiO4_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Zr.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.290) then
        call READ_DATEI('H2WO4.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('W.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.291) then
        call READ_DATEI('WCl.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.292) then
        call READ_DATEI('WF.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('F.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.293) then
        call READ_DATEI('WO2Cl2.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 2.D0
        stoich(4) = 2.D0
      elseif (specie.eq.294) then
        call READ_DATEI('W2O6.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 6.D0
      elseif (specie.eq.295) then
        call READ_DATEI('W3O8.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 3.D0
        stoich(3) = 8.D0
      elseif (specie.eq.296) then
        call READ_DATEI('W3O9.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 3.D0
        stoich(3) = 9.D0
      elseif (specie.eq.297) then
        call READ_DATEI('W4O12.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('W.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 4.D0
        stoich(3) = 12.D0
      elseif (specie.eq.298) then
        call READ_DATEI('Ni_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ni.txt'    ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.299) then
        call READ_DATEI('Ni_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ni.txt'    ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.300) then
        call READ_DATEI('Cr_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Cr.txt'    ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.301) then
        call READ_DATEI('Cr_l.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Cr.txt'    ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.302) then
        call READ_DATEI('CrN_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('CrN.txt'    ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.303) then
        call READ_DATEI('V2O3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('V.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 3.D0
      elseif (specie.eq.304) then
        call READ_DATEI('V2O4_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('V.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 4.D0
      elseif (specie.eq.305) then
        call READ_DATEI('V2O5_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('V.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 5.D0
      elseif (specie.eq.306) then
        call READ_DATEI('Ni3S2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ni.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 3.D0
        stoich(3) = 2.D0
      elseif (specie.eq.307) then
        call READ_DATEI('Ni3S2_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Ni.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 3.D0
        stoich(3) = 2.D0
      elseif (specie.eq.308) then
        call READ_DATEI('NH4Cl_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('N.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('Cl.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 4.D0
        stoich(4) = 1.D0
      elseif (specie.eq.309) then
        call READ_DATEI('P_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.310) then
        call READ_DATEI('P_l.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.311) then
        call READ_DATEI('P2O10_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 2.D0
        stoich(3) = 10.D0
      elseif (specie.eq.312) then
        call READ_DATEI('P4S3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 4.D0
        stoich(3) = 3.D0
      elseif (specie.eq.313) then
        call READ_DATEI('P4S3_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'     ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'     ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 4.D0
        stoich(3) = 3.D0
      elseif (specie.eq.314) then
        call READ_DATEI('VO.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('V.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      elseif (specie.eq.315) then
        call READ_DATEI('VO2.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('V.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 2.D0
      elseif (specie.eq.316) then
        call READ_DATEI('SO3.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('S.txt' ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt' ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.317) then
        call READ_DATEI('Zn_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Zn.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.318) then
        call READ_DATEI('Zn_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Zn.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.319) then
        call READ_DATEI('ZnSO4_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Zn.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('S.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.320) then
        call READ_DATEI('H3PO4_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'    ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('P.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 3.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.321) then
        call READ_DATEI('Mg3P2O8_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Mg.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('P.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 3.D0
        stoich(3) = 2.D0
        stoich(4) = 8.D0
      elseif (specie.eq.322) then
        call READ_DATEI('P3N5_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('N.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 3.D0
        stoich(3) = 5.D0
      elseif (specie.eq.323) then
        call READ_DATEI('P4O6.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 4.D0
        stoich(3) = 6.D0
      elseif (specie.eq.324) then
        call READ_DATEI('P4O10.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('P.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 4.D0
        stoich(3) = 10.D0
      elseif (specie.eq.325) then
        call READ_DATEI('AlF3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('F.txt'    ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 3.D0
      elseif (specie.eq.326) then
        call READ_DATEI('CaF2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('CaF2.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.327) then
        call READ_DATEI('KF_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('KF.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.328) then
        call READ_DATEI('NaF_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('NaF.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.329) then
        call READ_DATEI('H3PO4_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('H.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('P.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 3.D0
        stoich(3) = 1.D0
        stoich(4) = 4.D0
      elseif (specie.eq.330) then
        call READ_DATEI('FeF2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('FeF2.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.331) then
        call READ_DATEI('MgF2_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('MgF2.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.332) then
        call READ_DATEI('AlF6Na3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Al.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('F.txt'    ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('Na.txt'   ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 1.D0
        stoich(3) = 6.D0
        stoich(4) = 3.D0
      elseif (specie.eq.333) then
        call READ_DATEI('Li2SiO3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.334) then
        call READ_DATEI('Li2SiO3_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.335) then
        call READ_DATEI('Li2Si2O5_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 2.D0
        stoich(4) = 5.D0
      elseif (specie.eq.336) then
        call READ_DATEI('Li2Si2O5_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Si.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 2.D0
        stoich(4) = 5.D0
      elseif (specie.eq.337) then
        call READ_DATEI('Li2TiO3_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.338) then
        call READ_DATEI('Li2TiO3_l.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Li.txt'   ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('Ti.txt'   ,dG,T,Nmax,N,S,3) 
        call READ_DATEI('O.txt'    ,dG,T,Nmax,N,S,4) 
        Edzahl = 3
        stoich(2) = 2.D0
        stoich(3) = 1.D0
        stoich(4) = 3.D0
      elseif (specie.eq.339) then
        call READ_DATEI('Co_cr.txt',dG,T,Nmax,N,S,1) 
        call READ_DATEI('Co.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.340) then
        call READ_DATEI('Co_l.txt' ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Co.txt'   ,dG,T,Nmax,N,S,2) 
        Edzahl = 1
        stoich(2) = 1.D0
      elseif (specie.eq.341) then
        call READ_DATEI('CoO_cr.txt'  ,dG,T,Nmax,N,S,1) 
        call READ_DATEI('Co.txt'      ,dG,T,Nmax,N,S,2) 
        call READ_DATEI('O.txt'       ,dG,T,Nmax,N,S,3) 
        Edzahl = 2
        stoich(2) = 1.D0
        stoich(3) = 1.D0
      else
        write(*,*) 'Specie=',specie,' ???'
        stop
      endif
      write(*,'("stoichiometry:",99(I3))') 1,INT(stoich(2:Edzahl+1))
      write(*,*)
*
      write(*,*) 'T data from ... to ...'
      write(*,*) (T(i,1),i=1,1+Edzahl)
      write(*,*) (T(i,Nmax(i)),i=1,1+Edzahl)
      write(*,*)
      write(*,*) 'fit for which T-range ...?'
      read(*,*) T1,T2
      write(*,*)
      write(*,*) 'fit for what?'
      write(*,*)
      write(*,*) '0 = dG(T) without 1/T-term'
      write(*,*) '1 = dG(T)'
      write(*,*) '2 = kp(T) Gail-fit'
      write(*,*) '3 = -dG/RT Stock-fit'
      write(*,*) '4 = pvap(T)'
      read(*,*) mode
      write(*,*)
      Ndat = 0
      Tdatmin = 9.d+99
      Tdatmax = 0.d0
      sumst = 1
      do ed=2,Edzahl+1
        sumst = sumst - stoich(ed)    ! (1-Natom)
      enddo
      do i=1,Nmax(1)
        if ((T(1,i).ge.T1).and.(T(1,i).le.T2)) then
          xx = T(1,i)
          yy = dG(1,i)
          ok = .true.
          do ed=2,Edzahl+1
            gef = .false.
            do j=1,Nmax(ed)
              if ((.not.gef).and.(DABS(xx-T(ed,j)).lt.1.D-3)) then
                yy  = yy - stoich(ed)*dG(ed,j)
                gef = .true.
              endif
            enddo
            ok = (ok.and.gef)
          enddo
          ! xx is T [K]
          ! yy is \Delta_f G^1bar [kJ/mol]
          if (ok) then
            Ndat = Ndat + 1
            TTform(Ndat) = xx
            dGform(Ndat) = yy
            Tdatmin = MIN(Tdatmin,xx)
            Tdatmax = MAX(Tdatmax,xx)
            if (mode.eq.1) then
              x(Ndat) = xx
              y(Ndat) = yy * 1000.D0 * xx  
            else if (mode.eq.2) then
              x(Ndat) = 5040.D0/xx
              y(Ndat) = sumst*DLOG(bar) - yy*1000.D0/(Rgas*xx)
            else if (mode.eq.3) then
              x(Ndat) = xx
              y(Ndat) = -yy*1000.D0/(Rgas*xx)
            else if (mode.eq.4) then
              x(Ndat) = xx
              y(Ndat) = (DLOG(bar) + yy*1000.D0/(Rgas*xx))*xx
            else if (mode.eq.0) then
              x(Ndat) = xx
              y(Ndat) = yy * 1000.D0
            endif
            write(*,1010) xx,yy,x(Ndat),y(Ndat)
          endif
        endif
      enddo

      !--- linear data extrapolation for low T ---
      if (Tdatmin>T1+5.d0) then
        Ggrad = (dGform(3)-dGform(1))
     >         /(TTform(3)-TTform(1))
        write(*,*) "linear extrapolation dGform(T) ..."
        dN = (Tdatmin-T1)/20
        if (dN.eq.0) dN=1
        x(1+dN:Ndat+dN) = x(1:Ndat)
        y(1+dN:Ndat+dN) = y(1:Ndat)
        dGform(1+dN:Ndat+dN) = dGform(1:Ndat)
        TTform(1+dN:Ndat+dN) = TTform(1:Ndat)
        Ndat = Ndat+dN
        do i=1,dN
          xx = T1 + (Tdatmin-T1)*DBLE(i-1)/DBLE(dN)
          yy = dGform(1) + Ggrad*(xx-TTform(1))
          if (mode.eq.1) then
            x(i) = xx
            y(i) = yy * 1000.D0 * xx 
          else if (mode.eq.2) then
            x(i) = 5040.D0/xx
            y(i) = sumst*DLOG(bar) - yy*1000.D0/(Rgas*xx)
          else if (mode.eq.3) then
            x(i) = xx
            y(i) = -yy*1000.D0/(Rgas*xx)
          else if (mode.eq.4) then
            x(i) = xx
            y(i) = (DLOG(bar) + yy*1000.D0/(Rgas*xx))*xx
          else if (mode.eq.0) then
            x(i) = xx
            y(i) = yy * 1000.D0
          endif
          write(*,1010) xx,yy,x(i),y(i)
        enddo
      endif  

      !--- linear data extrapolation for high T ---
      if (Tdatmax<T2-5.d0) then
        Ggrad = (dGform(Ndat)-dGform(Ndat-1))
     >         /(TTform(Ndat)-TTform(Ndat-1))
        write(*,*) "linear extrapolation dGform(T) ..."
        dN = (T2-Tdatmax)/100
        if (dN.eq.0) dN=1
        do i=1,dN
          xx = Tdatmax + (T2-Tdatmax)*DBLE(i)/DBLE(dN)
          yy = dGform(Ndat) + Ggrad*(xx-TTform(Ndat))
          if (mode.eq.1) then
            x(Ndat+i) = xx
            y(Ndat+i) = yy * 1000.D0 * xx 
          else if (mode.eq.2) then
            x(Ndat+i) = 5040.D0/xx
            y(Ndat+i) = sumst*DLOG(bar) - yy*1000.D0/(Rgas*xx)
          else if (mode.eq.3) then
            x(Ndat+i) = xx
            y(Ndat+i) = -yy*1000.D0/(Rgas*xx)
          else if (mode.eq.4) then
            x(Ndat+i) = xx
            y(Ndat+i) = (DLOG(bar) + yy*1000.D0/(Rgas*xx))*xx
          else if (mode.eq.0) then
            x(Ndat+i) = xx
            y(Ndat+i) = yy * 1000.D0
          endif
          write(*,1010) xx,yy,x(Ndat+i),y(Ndat+i)
        enddo
        Ndat = Ndat+dN
      endif  

      write(*,*) Ndat,' data points'
 100  write(*,*)

      if (mode.eq.3) then
        grad = 1 
        Bfit(:) = 1.E-6
        Bfit(2) = -1.d0
        call POLYFIT(N,1,Ndat,x,(y-Bfit(2)*log(x))*x,grad,koeff)
        Bfit(1) = koeff(0)
        Bfit(3) = koeff(1)
        e   = 1.E-6
        ee1 = 0.1d0
        write(*,*) 'polynomial degree=? ...'
        read(*,*) grad
        NB=grad+3
        Bfit(NB+1:5) = 0.d0
        print'("coefficients before:",5(1pE18.10))',Bfit(1:NB)
        do it=1,100
          Nit = 0 
          call PARAM_LS(stock,NB,Nit,Ndat,dev,e,ee1,Bfit,x,y)
          print'("coefficients after :",5(1pE18.10))',Bfit(1:NB)
          print*,Nit,' iterations'
        enddo
        koeff(0:4) = Bfit(1:5)
        grad = 4
      else if ((mode.eq.4).or.(mode.eq.5)) then
        mode = 4 
        write(*,*) 'do you want ln pvap = A+B/(T+C) fit?  (y/n)'
        read(*,*) answer
        if (answer.eq.'y') then
          grad = 1 
          call POLYFIT(N,1,Ndat,x,y,grad,koeff)
          Afit(1) = koeff(1)
          Afit(2) = koeff(0) 
          Afit(3) = -5.d0
          e   = 1.E-6
          ee1 = 0.5d0
          do it=1,20
            print'("coefficients before:",5(1pE18.10))',Afit(1:3)
            call PARAM_LS(pvap,3,Nit,Ndat,dev,e,ee1,Afit,x,y/x)
            print'("coefficients after :",5(1pE18.10))',Afit(1:3)
            print*,Nit,' iterations'
          enddo  
          koeff(0:2) = Afit(1:3)
          grad = 2
          mode = 5
        endif
      endif
      if (mode<=2.or.mode==4) then
        write(*,*) 'polynomial degree=? ...'
        read(*,*) grad
        call POLYFIT(N,1,Ndat,x,y,grad,koeff)
      endif

      if (mode.eq.1) then
        write(*,*)
        write(*,*) 'dG(T)[J/mol] = Summe_i{ a_i*T^(i-1) }'
        write(*,*)
      else if (mode.eq.2) then
        write(*,*)
        write(*,*) 'kp(T)[cgs] = EXP( Summe_i{ a_i*(5040/T)^i } )'
        write(*,*)
      else if (mode.eq.3) then
        write(*,*)
        write(*,*) '-dG/RT = a0/T + a1*ln(T) + a2 + a3*T + a4*T^2'
        write(*,*)
      else if (mode.eq.4) then
        write(*,*)
        write(*,*) 'ln pvap(T)[dyn/cm2] = Summe_i{ a_i*T^(i-1) }'
        write(*,*)
      else if (mode.eq.5) then
        write(*,*)
        write(*,*) 'ln pvap(T)[dyn/cm2] = a_0 + a_1 / ( T + a_2 )'
        write(*,*)
      else if (mode.eq.0) then
        write(*,*)
        write(*,*) 'dG(T)[J/mol] = Summe_i{ a_i*T^i }'
        write(*,*)
      endif  
      do i=0,grad
        write(*,1000) i,koeff(i)
      enddo
      write(*,*)
      delta1 = 0.D0
      delta2 = 0.D0
      open(unit=12,file='fitout.dat',status='replace')
      do i=1,Ndat
        xx = x(i)
        yy = y(i)
        yfit = 0.D0
        do j=0,grad
          yfit = yfit + koeff(j)*xx**j
        enddo
        if (mode.eq.1) then
          yy   = yy/1000.D0/xx
          yfit = yfit/1000.D0/xx
        else if (mode.eq.2) then
        else if (mode.eq.3) then
          yy   = -yy*(Rgas*xx)/1000.0 
          yfit = -stock(NB,Bfit,xx)*(Rgas*xx)/1000.0 
        else if (mode.eq.4) then
          yy   = yy/xx
          yfit = yfit/xx
        else if (mode.eq.5) then
          yy   = yy/xx
          yfit = pvap(3,Afit,xx)
        else if (mode.eq.0) then
          yy   = yy/1000.D0
          yfit = yfit/1000.D0
        endif
        delta1 = MAX(delta1, ABS(yy-yfit))
        delta2 = delta2 + (yy-yfit)**2
        write(*,1020)  5040.d0/xx,xx,yy,yfit
        write(12,1020) xx,yy,yfit
      enddo
      close(12)
      open(unit=12,file='fitout2.dat',status='replace')
      do i=0,1000
        xx = EXP(LOG(50.d0)+LOG(10000.d0/50.d0)*i/1000.d0)
        if (mode.eq.2) xx=5040.d0/xx
        yfit = 0.d0
        do j=0,grad
          yfit = yfit + koeff(j)*xx**j
        enddo
        if (mode.eq.1) then
          yfit = yfit/1000.D0/xx
        else if (mode.eq.2) then
        else if (mode.eq.3) then
          yy   = -yy*(Rgas*xx)/1000.0 
          yfit = -stock(NB,Bfit,xx)*(Rgas*xx)/1000.0 
        else if (mode.eq.4) then
          yy   = yy/xx
          yfit = yfit/xx
        else if (mode.eq.5) then
          yy   = yy/xx
          yfit = pvap(3,Afit,xx)
        else if (mode.eq.0) then
          yfit = yfit/1000.D0
        endif
        write(12,1020) xx,yfit
      enddo
      close(12)
      delta2 = DSQRT(delta2/dble(Ndat))
      if (mode.eq.1) then
        write(*,1031) delta1,delta2
      else
        write(*,1032) delta1,delta2
      endif
      write(*,*)
      write(*,2000) (koeff(j),j=0,grad)
      goto 100
*
      STOP
 1000 format(' a_',i1,' = ',1pE12.5)
 1010 format(2(0pF9.2),99(1pE14.5))
 1020 format(99(0pF13.6))
 1031 format('Abweichung max/mean = ',1pE9.2,' kJ/mol',1pE9.2)
 1032 format('Abweichung max/mean = ',2(1pE9.2))
 2000 format(99(1pE13.5))
      end
*
*
* ------------------------------------------------
      SUBROUTINE READ_DATEI(datei,y,x,Nmax,N,S,Nr)
* ------------------------------------------------
      implicit none
      integer :: N,S,Nr,i,ok,Nmax(S)
      real*8  :: x(S,N),y(S,N),T,Cp,SS,HGoverT,HH,dH,dG,T1000
      real*8  :: rgas=8.31434D+0,bar=1.D+6
      character*(*) datei
      character*100 zeile

      Nmax(Nr) = 0
      write(*,*) 'Lese Datei ',datei,' ...'
      open(unit=12,file=datei,status='old')
      do
        read(12,'(A100)',end=999) zeile
        read(zeile,*,iostat=ok) T,Cp,SS,HGoverT,HH,dH,dG
        !write(*,*) trim(zeile),ok
        if (ok==0) then
          Nmax(Nr) = Nmax(Nr) + 1
          x(Nr,Nmax(Nr)) = T              ! [K]
          y(Nr,Nmax(Nr)) = dG             ! [kJ/mol]
          !T1000 = T/1000.d0
          !y(Nr,Nmax(Nr)) = dH - HGoverT*T1000 -rgas*T1000
        endif   
      enddo
 999  close(12)    
      end

* -------------------------------------------------
      SUBROUTINE READ_DATEI2(datei,y,x,Nmax,N,S,Nr)
* -------------------------------------------------
      implicit none
      integer N,S,Nr,i
      integer Nmax(S)
      real*8 x(S,N),y(S,N),T,dG,CTOF
      character*(*) datei
      character*100 zeile
      real*8 data(100)
*
      Nmax(Nr) = 0
      write(*,*) 'Lese Datei ',datei,' ...'
      open(unit=12,file=datei,status='old')
      do i=1,2
        read(12,1000) zeile
      enddo  
 10   continue
        read(12,*,end=999) data(1:8)
        T = data(1)
        dG = data(7)
        write(*,*) T,dG
        if ((T.gt.0.D0).and.(T.lt.6001.D0)) then
          Nmax(Nr) = Nmax(Nr) + 1
          x(Nr,Nmax(Nr)) = T
          y(Nr,Nmax(Nr)) = dG
        endif
      if (T.lt.5999.D0) goto 10
 999  continue
      close(12)    
      RETURN
 1000 format(a100)
      end

* ----------------------------------
      INTEGER FUNCTION CTOI(N1,N2,c)
* ----------------------------------
      implicit none
      integer N1,N2,i,wert,minus
      character c(N2)
      wert = 0
      minus = 1
      do i=N1,N2
        if (c(i).eq.'-') then
          minus = -1
        else if (c(i).ge.'0') then
	  wert = 10*wert + ICHAR(c(i)) - ICHAR('0')
        else if (c(i).ne.' ') then
          write(*,*) 'komisches Zeichen in ctoi',c(i)
          STOP
	endif
      enddo
      CTOI = minus*wert
      end

* ---------------------------------
      REAL*8 FUNCTION CTOF(N1,N2,c)
* ---------------------------------      
      implicit none
      integer N1,N2,j,j2,CTOI,INSTR
      character c(N2)
      integer vor,nach
      real*8  pot
      j = INSTR(N1,N2,c,'.')
      if (j.gt.0) then
        vor  = CTOI(N1,j-1,c)
        nach = CTOI(j+1,N2,c)
        j2   = instr(j+1,N2,c,' ')
        if (j2.eq.0) j2=N2+1
        pot = 10.D0**DBLE(j2-j-1)
        if ((pot.le.1.D0).and.(nach.gt.0)) then
          write(*,*) '*** ',c,pot
	  STOP
	endif
        CTOF = DBLE(vor) + DSIGN(DBLE(nach)/pot,DBLE(vor))
      else
        CTOF = DBLE(CTOI(N1,N2,c))
      endif	
      end

* ----------------------------------------
      INTEGER FUNCTION INSTR(N1,N2,c,such)
* ----------------------------------------      
      implicit none
      integer N1,N2,i
      character c(N2),such
      INSTR=0
      do i=N1,N2
        if (c(i).eq.such) then
          INSTR = i
	  goto 10
	endif
      enddo
   10 continue    
      end
*
*
***********************************************************************
      SUBROUTINE POLYFIT(N,N1,N2,x,y,grad,koeff)
***********************************************************************
      implicit none
      integer N,N1,N2,grad,i,j,gmax
      parameter (gmax=20)
      real*8 x(N),y(N),koeff(grad+1)
      real*8 sumxi(0:2*gmax),sumxiyi(0:gmax)
      real*8 AA(gmax,gmax),bb(gmax),xx(gmax),ly(1000),dx(1000)
      real*8 A0(gmax,gmax),b0(gmax),sum
*
      do j=0,2*grad
        sumxi(j) = 0.D0
        do i=N1,N2
          sumxi(j) = sumxi(j) + x(i)**j
	enddo
      enddo	
      do j=0,grad
        sumxiyi(j) = 0.D0
        do i=N1,N2
          sumxiyi(j) = sumxiyi(j) + y(i)*x(i)**j
	enddo  
      enddo
      do i=0,grad
        do j=0,grad
          AA(i+1,j+1) = sumxi(i+j)
        enddo
        bb(i+1) = sumxiyi(i)
      enddo
      call GAUSS2(gmax,AA,xx,bb,grad+1)
      do i=1,grad+1
        koeff(i) = xx(i)
      enddo
*      
      RETURN
 1000 format(99(1pD23.15))
      end
*
*      
**********************************************************************
      SUBROUTINE GAUSS2 (dim,a,x,b,N)
**********************************************************************
*****                                                            *****
*****   wie GAUSS - aber mit Uebergabe von dim*dim-Matrix und    *****
*****   dim-Vektoren, von denen nur jeweils die ersten N Werte   *****
*****   interessieren.                                           *****
*****                                                            *****
**********************************************************************
*
      implicit none
      integer N,dim,i,j,k,kmax
      real*8  a(dim,dim),x(dim),b(dim),c,amax
*
      do 500 i=1,N-1
*       ------------------------------------------
*       ***  MAX-Zeilentausch der i-ten Zeile  ***      
*       ------------------------------------------
        kmax = i
        amax = DABS(a(i,i))
        do 200 k=i+1,N
          if ( DABS(a(k,i)) .gt. amax ) then
            kmax = k
            amax = DABS(a(k,i))
          endif
  200   continue
        if (kmax.ne.i) then
          do 210 j=1,N
            c         = a(i,j)
            a(i,j)    = a(kmax,j)
            a(kmax,j) = c 
  210     continue
          c       = b(i)
          b(i)    = b(kmax)
          b(kmax) = c 
        endif
*
*       ---------------------------------
*       ***  bringe auf Dreiecksform  ***
*       ---------------------------------
        do 310 k=i+1,N
          c = a(k,i) / a(i,i)
          a(k,i) = 0.D0
          do 300 j=i+1,N
            a(k,j) = a(k,j) - c * a(i,j)
  300     continue        
          b(k) = b(k) - c * b(i)
  310   continue
*
  500 continue
*
*     --------------------------
*     ***  loese nach x auf  ***
*     --------------------------
      do 610 i=N,1,-1
        c = 0.D0
        if (i.lt.N) then
          do 600 j=i+1,N
            c = c + a(i,j) * x(j)
  600     continue
        endif
        x(i) = (b(i) - c) / a(i,i)
  610 continue
*
      RETURN
      end

***************************************************************
      REAL*8 FUNCTION stock(Npara,Bfit,T)
***************************************************************
      implicit none
      integer,intent(in) :: Npara
      real*8,intent(in)  :: Bfit(Npara),T
      stock = Bfit(1)/T + Bfit(2)*LOG(T) + Bfit(3)
     >      + Bfit(4)*T + Bfit(5)*T**2
      end

***************************************************************
      REAL*8 FUNCTION pvap(Npara,Afit,T)
***************************************************************
      implicit none
      integer,intent(in) :: Npara
      real*8,intent(in)  :: Afit(Npara),T
      pvap = Afit(1)+Afit(2)/(T+Afit(3))
      end

***************************************************************
      SUBROUTINE S200(func,Npara,l2,Ndat,d,A,X,Y) 
***************************************************************
      implicit none
      integer,parameter :: SIZE=200
      integer,intent(in) :: Npara,Ndat
      real*8 :: d,l2,A(Npara),X(SIZE),Y(SIZE) 
      real*8 :: xx,yy
      real*8,external :: func
      integer :: j
      l2 = 0.d0
      do j = 1,Ndat
         xx = X(j)
         yy = func(Npara,A,xx)
         l2 = l2 + (Y(j)-yy)**2
      enddo
      d = dsqrt(l2/dfloat(Ndat-Npara))
      end

!***************************************************************
!* Parametric least squares curve fit subroutine. This program *
!* least squares fits a function to a set of data values by    *
!* successively reducing the variance. Convergence depends on  *
!* the initial values and is not assured.                      *
!* n pairs of data values, X(i), Y(i), are given. There are l  *
!* parameters, A(j), to be optimized across.                   *
!* Required are initial values for the A(l) and e. Another     *
!* important parameter which affects stability is e1, which is *
!* initially converted to e1(l) for the first intervals.       *
!* The parameters are multiplied by (1 - e1(i)) on each pass.  *
!***************************************************************
      SUBROUTINE PARAM_LS(func,l,m,n,d,e,ee1,A,X,Y)  
      implicit none
      integer,parameter   :: SIZE=200
      integer,intent(in)  :: l,n
      integer,intent(out) :: m
      real*8,intent(in)   :: e,ee1,X(SIZE),Y(SIZE)
      real*8,intent(out)  :: d,A(l)
      real*8,external :: func
      real*8  :: a0,l1,l2,m0,m1
      real*8  :: E1(l)
      integer :: i
      do i = 1,l
         E1(i) = ee1
      end do	
      !Set up test residual
      l1 = 1.d6
      !Make sweep through all parameters
 50   do i = 1,l
         a0 = A(i)
         !Get value of residual
         A(i) = a0
 100     call S200(func,l,l2,n,d,A,X,Y)
         !Store result in m0
         m0 = l2
         !Repeat for m1
         A(i) = a0 * (1.d0 - E1(i))
         call S200(func,l,l2,n,d,A,X,Y)
         m1 = l2
         !Change interval size if called for
         !If variance was increased, halve E1(i) 
         if (m1 > m0)  then
            E1(i) = -E1(i) / 2.d0
         end if
         !If variance was reduced, increase step size by increasing E1(i)
         if (m1 < m0)  then
            E1(i) = 1.2d0 * E1(i)
         end if 
         !If variance was increased, try to reduce it
         if (m1 > m0) then  
            A(i) = a0
         end if 
         if (m1 > m0) then
            goto 100
         end if
      enddo 
      !End of a complete pass
      !Test for convergence 
      m = m + 1
      if (l2.eq.0.d0) then
         return
      end if 
      if (dabs((l1-l2)/l2)<e)  then
         return
      end if 
      !If this point is reached, another pass is called for
      l1 = l2
      goto 50
      end
