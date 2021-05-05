{
    TFile* file_0 = new TFile("templates_521.root", "READ");
    TFile* file_1 = new TFile("templates_522.root", "READ");
    TFile* out_file= new TFile("templates.root", "RECREATE");

    for ( int i = 0; i < 712; i++)
    {
        string name_0=Form("Template_Ch%d.root",i);
        TH1D* hist_0 = (TH1D*)file_0->Get(name_0.c_str());

        if (hist_0->GetMean() > 0.0)
        {
            out_file-cd();
            hist_0.Write();
        }else
        {
            string name_1=Form("Template_Ch%d.root",i);
            TH1D* hist_1 = (TH1D*)file_1->Get(name_1.c_str());
            out_file-cd();
            hist_1.Write();

            delete hist_1;
        }

        delete hist_0;
    }
    file_0->Close();
    file_1->Close();
    out_file->Close();
}