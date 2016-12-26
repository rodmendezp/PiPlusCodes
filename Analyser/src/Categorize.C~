TString TIdentificator::GetCategorization(Int_t k)
{
    Int_t number_dc = fCT->GetNRows("DCPB");
    Int_t number_cc = fCT->GetNRows("CCPB");
    Int_t number_sc = fCT->GetNRows("SCPB");
    Int_t number_ec = fCT->GetNRows("ECPB");

    TString partId;

    partId = "not recognized";

    if (number_dc != 0) {
        if (k == 0 &&
                    Status(0) > 0 && Status(0) < 100 &&
                    Charge(0) == -1 &&
                    number_cc != 0 && number_ec != 0 && number_sc != 0 &&
                    StatCC(0) > 0 && StatSC(0) > 0 &&
                    StatDC(0) > 0 && StatEC(0) > 0 &&
                    DCStatus(0) > 0 && SCStatus(0) == 33 &&
                    Nphe(0) > 25 &&
                    Etot(0) / 0.27 / 1.15 + 0.4 > Momentum(0) &&
                    Etot(0) / 0.27 / 1.15 - 0.2 < Momentum(0) &&
                    Ein(0) + Eout(0) > 0.8 * 0.27 * Momentum(0) &&
                    Ein(0) + Eout(0) < 1.2 * 0.27 * Momentum(0) &&
                    Eout(0) != 0 && FidCheckCut() == 1)
            partId = "electron";


        if (k > 0) {
            if (Charge(k)==0 && Betta(k)>0.95 && ECStatus(k)>0)
                partId = "photon";

            //positive particles
            if (Charge(k) == 1 &&
                        Status(k) > 0 && Status(k) < 100 &&
                        StatDC(k) > 0 && DCStatus(k) > 0) {
                if (Momentum(k)>=2.7) {
                    if (number_cc != 0 && StatCC(k) > 0 &&
                                Nphe(k) > 25 && Chi2CC(k) < 5 / 57.3)
                        partId = "high energy pion +";
                }

                if (Momentum(k) < 2.7) {
                    if (number_sc != 0 && StatSC(k) > 0 &&
                              	((Momentum(k) > 0 && Momentum(k) <= 0.25 && TimeCorr4(0.139,k) >= -1.45 && TimeCorr4(0.139,k) <= 1.05) ||
                                (Momentum(k) > 0.25 && Momentum(k) <= 0.5 && TimeCorr4(0.139,k) >= -1.44 && TimeCorr4(0.139,k) <= 1.05) || 
				(Momentum(k) > 0.5 && Momentum(k) <= 0.75 && TimeCorr4(0.139,k) >= -1.44 && TimeCorr4(0.139,k) <= 1.05) ||
				(Momentum(k) > 0.75 && Momentum(k) <= 1 && TimeCorr4(0.139,k) >= -1.4 && TimeCorr4(0.139,k) <= 1.05) ||
				(Momentum(k) > 1 && Momentum(k) <= 1.25 && TimeCorr4(0.139,k) >= -1.35 && TimeCorr4(0.139,k) <= 1.03) ||
				(Momentum(k) > 1.25 && Momentum(k) <= 1.5 && TimeCorr4(0.139,k) >= -1.35 && TimeCorr4(0.139,k) <= 0.95) ||
				(Momentum(k) > 1.5 && Momentum(k) <= 1.75 && TimeCorr4(0.139,k) >= -1.35 && TimeCorr4(0.139,k) <= 0.87) ||
				(Momentum(k) > 1.75 && Momentum(k) <= 2 && TimeCorr4(0.139,k) >= -1.25 && TimeCorr4(0.139,k) <= 0.68) ||
				(Momentum(k) > 2 && Momentum(k) <= 2.25 && TimeCorr4(0.139,k) >= -0.95 && TimeCorr4(0.139,k) <= 0.65) ||
				(Momentum(k) > 2.25 && Momentum(k) <= 2.5 && TimeCorr4(0.139,k) >= -1.05 && TimeCorr4(0.139,k) <= 0.61
				&& Mass2(k) < 0.5) ||
				(Momentum(k) > 2.5 && Momentum(k) < 2.7 && TimeCorr4(0.139,k) >= -1.05 && TimeCorr4(0.139,k) <= 0.61
				&& Mass2(k) < 0.4)))
                        partId = "low energy pion +";
                }

                if (Momentum(k) < 2.) {
                    if (number_sc != 0 && StatSC(k) > 0 &&
                                ((Momentum(k) >= 1. &&
                                    TimeCorr4(0.938,k) >= -0.69 &&
                                    TimeCorr4(0.938,k) <= 1.38) ||
                                (Momentum(k) < 1. &&
                                    TimeCorr4(0.938,k) >= -3.78 &&
                                    TimeCorr4(0.938,k) <= 6.75)))
                        partId = "low energy proton";
                }

                if (Charge(k) == 1 && number_cc != 0 &&
                            number_ec != 0 && number_sc != 0 &&
                            StatCC(k) > 0 && StatSC(k) > 0 &&
                            StatDC(k) > 0 && StatEC(k) > 0 &&
                            DCStatus(k) > 0 && Nphe(k) > 25 &&
                            Etot(k) / 0.27 / 1.15 + 0.4 > Momentum(k) &&
                            Etot(k) / 0.27 / 1.15 - 0.2 < Momentum(k) &&
                            Ein(k) + Eout(k) > 0.8 * 0.27 * Momentum(k) &&
                            Ein(k) + Eout(k) < 1.2 * 0.27 * Momentum(k))
                    partId = "positron";
            }
        }
    }

    return partId;
}

TString TIdentificator::GetCategorizationOld(Int_t k){
    Int_t number_dc = fCT->GetNRows("DCPB");
    Int_t number_cc = fCT->GetNRows("CCPB");
    Int_t number_sc = fCT->GetNRows("SCPB");
    Int_t number_ec = fCT->GetNRows("ECPB");

    TString partId;

    partId = "not recognized";

    if (number_dc != 0) {
        if (k == 0 &&
                    Status(0) > 0 && Status(0) < 100 &&
                    Charge(0) == -1 &&
                    number_cc != 0 && number_ec != 0 && number_sc != 0 &&
                    StatCC(0) > 0 && StatSC(0) > 0 &&
                    StatDC(0) > 0 && StatEC(0) > 0 &&
                    DCStatus(0) > 0 && SCStatus(0) == 33 &&
                    Nphe(0) > 25 &&
                    Etot(0) / 0.27 + 0.4 > Momentum(0) &&
                    Etot(0) / 0.27 - 0.4 < Momentum(0) &&
                    Ein(0) + Eout(0) > 0.8 * 0.27 * Momentum(0) &&
                    Ein(0) + Eout(0) < 1.2 * 0.27 * Momentum(0) &&
                    FidCheckCut() == 1)
            partId = "electron";


        if (k > 0) {
            if (Charge(k)==0 && Betta(k)>0.95 && ECStatus(k)>0)
                partId = "photon";

            //positive particles
            if (Charge(k) == 1 &&
                        Status(k) > 0 && Status(k) < 100 &&
                        StatDC(k) > 0 && DCStatus(k) > 0) {
                if (Momentum(k)>=2.7) {
                    if (number_cc != 0 && StatCC(k) > 0 &&
                                Nphe(k) > 25 && Chi2CC(k) < 5 / 57.3)
                        partId = "high energy pion +";
                }

                if (Momentum(k) < 2.7) {
                    if (number_sc != 0 && StatSC(k) > 0 &&
                                ((Momentum(k) < 1 &&
                                    TimeCorr4(0.139,k) >= -1.46 &&
                                    TimeCorr4(0.139,k) <= 0.15) ||
                                (Momentum(k) >=1 &&
                                    TimeCorr4(0.139,k) >= -1.38 &&
                                    TimeCorr4(0.139,k) <= 0.53)))
                        partId = "low energy pion +";
                }

                if (Momentum(k) < 2.) {
                    if (number_sc != 0 && StatSC(k) > 0 &&
                                ((Momentum(k) >= 1. &&
                                    TimeCorr4(0.938,k) >= -0.69 &&
                                    TimeCorr4(0.938,k) <= 1.38) ||
                                (Momentum(k) < 1. &&
                                    TimeCorr4(0.938,k) >= -3.78 &&
                                    TimeCorr4(0.938,k) <= 6.75)))
                        partId = "low energy proton";
                }

                if (Charge(k) == 1 && number_cc != 0 &&
                            number_ec != 0 && number_sc != 0 &&
                            StatCC(k) > 0 && StatSC(k) > 0 &&
                            StatDC(k) > 0 && StatEC(k) > 0 &&
                            DCStatus(k) > 0 && Nphe(k) > 25 &&
                            Etot(k) / 0.27 + 0.4 > Momentum(k) &&
                            Etot(k) / 0.27 - 0.4 < Momentum(k) &&
                            Ein(k) + Eout(k) > 0.8 * 0.27 * Momentum(k) &&
                            Ein(k) + Eout(k) < 1.2 * 0.27 * Momentum(k))
                    partId = "positron";
            }
        }
    }
    return partId;
}



TString* TIdentificator::GetCategorization()
{
    Int_t number = fCT->GetNRows("EVNT");

    if (fPartIds != 0) delete [] fPartIds;
    fPartIds = new TString[number];

    if (number != 0) {
        for (Int_t i = 0; i < number; i++)
            fPartIds[i] = GetCategorization(i);
    }

    return fPartIds;
}



void TIdentificator::PrintCategorization()
{
    Int_t number = fCT->GetNRows("EVNT");
    fPartIds = GetCategorization();

    if (fPartIds->CompareTo("electron") == 0) {
        for (Int_t i = 0; i < number; i++)
            cout << *(fPartIds+i) << endl;
        cout << endl;
    }
}



void TIdentificator::PrintCategorization(TString* partIds)
{
    Int_t number = fCT->GetNRows("EVNT");

    if (partIds->CompareTo("electron") == 0) {
        for (Int_t i = 0; i < number; i++)
            cout << *(partIds+i) << endl;
        cout << endl;
    }
}

