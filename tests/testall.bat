python2.5 testcase.py -q "-mrx compliance" "-eSAS10mPRN7_4064.csv" "-R0.001" "-oSAS10mPRN7_4064"
python2.5 testcase.py -q "-mrx compliance" "-eSAS10mPRN7_6223.csv" "-R0.001" "-oSAS10mPRN7_6223"
python2.5 testcase.py -q "-mrx compliance" "-e4012_800mVPRBS7.csv" "-R0.001" "-o4012_800mVPRBS7"
python2.5 testcase.py -q "-mrx compliance" "-e4017_800mVPRBS7.csv" "-R0.001" "-o4017_800mVPRBS7"
python2.5 testcase.py -q "-mrx compliance" "-eSAS10m800mVPRBS7.csv" "-R0.001" "-oSAS10m800mVPRBS7"

python2.5 testcase.py -q "-mrx compliance" "-e0m-prbs10-12_5ps-ag.csv" "-o0m-prbs10-12_5ps-ag" "-R0.001" "-z5.0"
python2.5 testcase.py -q "-mrx compliance" "-e6m-prbs10-12_5ps-ag.csv" "-o6m-prbs10-12_5ps-ag" "-R0.001"
python2.5 testcase.py -q "-mrx compliance" "-e10m-prbs10-12_5ps-ag.csv" "-o10m-prbs10-12_5ps-ag" "-R0.001"
python2.5 testcase.py -q "-mrx compliance" "-eipass_10m_1vlaunch_prbs7_2.csv" "-oipass_10m_1vlaunch_prbs7_2" "-R0.001"
python2.5 testcase.py -q "-mrx compliance" "-eipass_10m_1vlaunch_prbs7_ffe1250_250.csv" "-oipass_10m_1vlaunch_prbs7_ffe1250_250" "-R0.001"

python2.5 testcase.py -q "-mrx compliance" "-eRxi_10m_prbs10_13.2sample165_march08.txt" "-x12.5e-12" "-R0.001" "-oRxi_10m_prbs10_13.2sample165_march08"
python2.5 testcase.py -q "-mtx compliance" "-eTxi_10m_prbs10_13.2sample165_march08.txt" "-x12.5e-12" "-z5.0" "-ssixmeter_A5A6B5B6.s4p" "-R0.001" "-oTxi_10m_prbs10_13.2sample165_march08"
python2.5 testcase.py -q "-shalfmeter_A5A6B5B6.s4p" "-ohalfmeter_A5A6B5B6"
python2.5 testcase.py -q "-ssixmeter_A5A6B5B6.s4p" "-osixmeter_A5A6B5B6"
python2.5 testcase.py -q "-sHP01_BtoB_3Connector.s4p" "-oHP01_BtoB_3Connector"
python2.5 testcase.py -q "-sSAS2_transmittertestload.s4p" "-oSAS2_transmittertestload"

