from lfsr import LFSR
from msalgo import MS
import random
import time
import pymysql
import func_timeout
# import timeout_decorator

def combine32_30(x1, x2, x3, x4, x5):
    return (x1 & x2) | (x1 & x3) | (x1 & x4) | ((x2 ^ 1) & (x3 ^ 1) & (x4 ^ 1) & x5)

def combine32_28(x1, x2, x3, x4, x5):
    return (x1 & x2) | (x1 & x3) | (x1 & x4) | ((x2 ^ 1) & (x3 ^ 1) & x4)

def combine32_26(x1, x2, x3, x4, x5):
    return (x1 & x2) | (x1 & x3) | ((x2 ^ 1) & x4 & x5) | ((x3 ^ 1) & x4 & x5)

def combine32_24(x1, x2, x3, x4, x5):
    return (x1 & x2) | (x1 & x3) | (x2 & (x3 ^ 1))

def combine32_22(x1, x2, x3, x4, x5):
    return (x1 & x2) | (x2 & (x3 ^ 1)) | (x1 & x3 & x5) | ((x3 ^ 1) & x4 & x5)

def combine32_20(x1, x2, x3, x4, x5):
    return (x1 & x2) | (x2 & (x3 ^ 1)) | ((x3 ^ 1) & x5)

def combine32_18(x1, x2, x3, x4, x5):
    return (x1 & x2) | (x2 & (x3 ^ 1)) | (x2 & x5) | (x2 & x4 & (x5 ^ 1))

db = pymysql.connect(host = 'localhost', user = 'root', password = 'root', db = 'msalgo', charset = 'utf8mb4')

data = {
    15 : 
        {
            2 : [16385,24576,16392,17408,16448,16512],
            4 : [16395,29696,16406,23040,16410,22016,16451,28800,16457,25728,16466,21120,16481,24960,16522,21568,16530,21056,16560,17216,16578,20672,16592,17088,16649,25632,16674,20768,16708,18592,16908,19472,16913,25104,16920,17936,16961,24720,16976,17040,17156,18480,17160,17456,17444,18696,17473,24712,17538,20552,17540,18504,18437,26628,18444,19460,18450,20996,18465,24836,18500,18564,18562,20548,18689,24612,20485,26626,20497,25090,20498,20994,20513,24834,20609,24642,20993,24594,24579,28673,24585,25601],
            6 : [16431,32000,16443,30464,16487,31104,16494,23936,16499,29568,16506,22400,16555,30016,16558,23872,16566,23360,16595,29376,16604,20160,16613,27072,16614,22976,16669,28192,16711,30880,16717,27808,16741,27040,16748,19872,16779,29792,16781,27744,16789,27232,16793,26208,16803,29024,16806,22880,16817,25440,16820,19296,16824,18272,16837,26848,16844,19680,16866,20960,16872,17888,16919,31248,16923,30224,16947,29456,16950,23312,16956,20240,16971,29840,16986,22160,17038,23632,17043,29264,17049,26192,17059,29008,17061,26960,17080,18256,17094,22736,17105,25296,17112,18128,17122,20944,17124,18896,17166,23600,17173,27184,17194,21808,17196,19760,17202,21296,17236,19120,17249,25008,17285,26736,17286,22640,17290,21616,17300,19056,17304,18032,17346,20720,17348,18672,17435,30216,17454,23816,17468,20232,17479,30856,17498,22152,17507,29064,17516,19848,17547,29768,17555,29256,17586,21320,17612,19656,17640,17864,17692,20008,17706,21800,17713,25384,17714,21288,17720,18216,17731,28840,17733,26792,17752,18088,17797,26728,17798,22632,17825,24936,17939,29208,17955,28952,17962,21784,17969,25368,17972,19224,17987,28824,17993,25752,18001,25240,18004,19096,18020,18840,18051,28760,18053,26712,18066,21080,18081,24920,18114,20696,18116,18648,18181,26680,18447,31748,18459,30212,18462,24068,18475,29956,18486,23300,18490,22276,18509,27780,18518,23172,18571,29764,18573,27716,18574,23620,18598,22852,18636,19652,18641,25284,18644,19140,18701,27684,18707,29220,18713,26148,18726,22820,18729,25892,18755,28836,18762,21668,18769,25252,18788,18852,18825,25700,18850,20836,18881,24804,18951,30740,18965,27156,18970,22036,18981,26900,18985,25876,19017,25748,19020,19604,19026,21140,19028,19092,19075,28756,19082,21588,19106,20820,19209,25652,19210,21556,19233,24884,19491,28940,19530,21644,19532,19596,19554,20876,19589,26700,19602,21068,19617,24908,19618,20812,19730,21036,19978,21532,20034,20636,20097,24668,20225,24636,20503,31234,20525,27906,20526,23810,20531,29442,20533,27394,20534,23298,20570,22146,20582,22914,20615,30786,20622,23618,20627,29250,20633,26178,20634,22082,20643,28994,20706,20930,20743,30754,20771,28962,20785,25378,20834,20898,20867,28770,20869,26722,20870,22626,20881,25186,21033,25874,21061,26770,21074,21138,21089,24978,21137,25170,21185,24786,21257,25650,21525,27146,21529,26122,21541,26890,21542,22794,21546,21770,21577,25738,21585,25226,21637,26698,21641,25674,21766,22570,21793,24874,21889,24682,22019,28698,22022,22554,22535,30726,22541,27654,22550,23046,22601,25734,22625,24966,22659,28742,22661,26694,22789,26662,22793,25638,22817,24870,23073,24854,23681,24654,24591,31745,24633,26369,24651,29825,24659,29313,24677,27009,24689,25473,24729,26177,24777,25793,24801,25025,24843,29729,25113,26129,25123,28945,25125,26897,25161,25745,25221,26705,25613,27657,26691,28805,26693,26757,26883,28709,28679,30723,28707,28931],
            8 : [16623,32192,16638,24512,16703,32544,16815,32096,16855,31456,16862,24288,17007,32144,17019,30608,17022,24464,17071,32080,17085,28496,17133,28112,17142,23504,17145,26576,17214,24368,17239,31408,17243,30384,17261,28080,17333,27504,17351,30960,17365,27376,17366,23280,17369,26352,17400,18416,17623,31432,17639,31176,17646,24008,17654,23496,17725,28456,17743,31912,17751,31400,17774,23976,17821,28264,17846,23400,17849,26472,17870,23784,17893,27112,17898,21992,17906,21480,17912,18408,18013,28312,18093,27992,18125,27864,18149,27096,18150,23000,18162,21464,18164,19416,18191,31800,18205,28216,18219,30008,18233,26424,18234,22328,18247,30904,18251,29880,18253,27832,18281,26040,18290,21432,18311,30840,18323,29304,18325,27256,18332,20088,18354,21368,18386,21240,18404,18936,18621,28484,18654,24260,18663,31172,18719,32292,18749,28452,18779,30372,18798,23972,18810,22436,18812,20388,18859,30052,18873,26468,18894,23780,18906,22244,18917,27108,19006,24340,19037,28308,19051,30100,19053,28052,19059,29588,19066,22420,19102,24148,19125,27476,19149,27860,19158,23252,19188,19412,19227,30260,19230,24116,19251,29492,19260,20276,19275,29876,19286,23220,19301,27060,19349,27252,19354,22132,19366,22900,19377,25460,19402,21748,19409,25332,19410,21236,19543,31372,19549,28300,19563,30092,19566,23948,19607,31308,19611,30284,19638,23372,19661,27852,19676,20172,19685,27084,19690,21964,19697,25548,19751,31020,19755,29996,19757,27948,19758,23852,19766,23340,19772,20268,19787,29868,19790,23724,19797,27308,19804,20140,19811,29100,19825,25516,19853,27756,19865,26220,19881,25964,19909,26860,19913,25836,19991,31260,19997,28188,20014,23836,20019,29468,20026,22300,20028,20252,20074,21916,20082,21404,20138,21852,20145,25436,20165,26844,20178,21212,20235,29756,20259,28988,20266,21820,20273,25404,20305,25276,20357,26748,20361,25724,20369,25212,20591,32130,20605,28546,20663,31554,20669,28482,20687,31938,20711,31170,20791,31522,20843,30114,20845,28066,20854,23458,20917,27490,20921,26466,20935,30946,20939,29922,20953,26338,20966,23010,20970,21986,21079,31378,21086,24210,21101,28050,21110,23442,21159,31058,21203,29394,21222,22994,21271,31282,21275,30258,21305,26418,21325,27826,21326,23730,21347,29106,21389,27762,21397,27250,21446,22770,21566,24330,21583,31882,21595,30346,21607,31114,21614,23946,21690,22346,21715,29386,21718,23242,21783,31274,21799,31018,21838,23722,21846,23210,21865,26026,21895,30826,21902,23658,21914,22122,21923,29034,21937,25450,21957,26858,22043,30234,22061,27930,22070,23322,22073,26394,22087,30874,22091,29850,22115,29082,22157,27738,22225,25306,22286,23610,22294,23098,22321,25402,22342,22714,22353,25274,22409,25722,22433,24954,22559,32262,22587,30470,22615,31366,22649,26502,22686,24134,22699,30022,22702,23878,22727,30918,22731,29894,22734,23750,22755,29126,22823,31014,22861,27814,22923,29798,22931,29286,22985,25830,23083,29974,23094,23318,23123,29334,23182,23638,23203,29014,23209,25942,23265,25046,23307,29750,23309,27702,23315,29238,23321,26166,23333,26934,23365,26806,23489,24822,23591,30990,23630,23694,23657,25998,23699,29262,23701,27214,23715,29006,23777,25038,23875,28846,23881,25774,23953,25198,23969,24942,24089,26142,24133,26782,24161,24990,24201,25694,24257,24798,24337,25150,24385,24766,24701,28545,24735,32321,24751,32065,24763,30529,24765,28481,24813,28097,24819,29633,25003,30049,25005,28001,25017,26465,25037,27873,25059,29153,25073,25569,25135,32017,25143,31505,25167,31889,25179,30353,25205,27537,25239,31313,25261,27985,25293,27857,25305,26321,25397,27441,25415,30897,25419,29873,25429,27313,25443,29105,25513,25969,25647,32009,25743,31817,25755,30281,25757,28233,25773,27977,25827,29129,25829,27081,25879,31273,25907,29481,25927,30889,25995,29801,25997,27753,26005,27241,26051,28905,26189,27801,26195,29337,26201,26265,26261,27225,26307,28889,26387,29241,26655,32261,26683,30469,26727,31109,26731,30085,26775,31301,26797,27973,26829,27845,26907,30245,26923,29989,26933,27429,26951,30885,27019,29797,27027,29285,27165,28181,27189,27413,27211,29845,27221,27285,27271,30805,27283,29269,27299,29013,27411,29237,27691,29965,27719,30861,27725,27789,27731,29325,27811,29005,27915,29741,27923,29229,28419,28733,28719,32003,28763,30339,28779,30083,28839,31043,28851,29507,28967,31011,28979,29475,29063,30819,29199,31763,29259,29843,29323,29779,29719,31243,29831,30795,30735,31751,30743,31239],
            10 : [17279,32688,17375,32496,17391,32240,17663,32712,18111,32600,18143,32472,18167,31704,18271,32440,18301,28600,18363,30584,18366,24440,18879,32612,18927,32228,18935,31716,18941,28644,19071,32660,19183,32212,19311,32180,19325,28596,19407,31988,19437,28148,19438,24052,19449,26612,19647,32588,19679,32460,19775,32556,19837,28588,19887,32108,19895,31596,19899,30572,19931,30444,19934,24300,19957,27628,20063,32412,20091,30620,20175,31964,20190,24284,20303,31932,20342,23484,20345,26556,20375,31356,20391,31100,20403,29564,20429,27900,20438,23292,20441,26364,20442,22268,20454,23036,21374,24498,21423,32114,21469,28402,21470,24306,21479,31218,21486,24050,21491,29682,21751,31690,21755,30666,21981,28394,21995,30186,21997,28138,22009,26602,22127,32154,22142,24474,22199,31578,22231,31450,22254,24026,22259,29658,22261,27610,22303,32314,22319,32058,22331,30522,22365,28346,22379,30138,22382,23994,22387,29626,22394,22458,22415,31866,22423,31354,22429,28282,22451,29562,22471,30970,22485,27386,22505,26106,22513,25594,22767,32198,22847,32550,22895,32166,22903,31654,22907,30630,22973,28518,23003,30438,23019,30182,23021,28134,23022,24038,23103,32534,23135,32406,23159,31638,23163,30614,23199,32342,23227,30550,23229,28502,23261,28374,23275,30166,23283,29654,23285,27606,23357,28470,23383,31414,23411,29622,23413,27574,23439,31862,23467,30070,23477,27510,23502,23798,23513,26358,23647,32398,23677,28558,23739,30542,23767,31438,23773,28366,23801,26574,23899,30382,23915,30126,23966,24174,23981,28014,23987,29550,24037,27118,24125,28446,24143,31902,24215,31326,24243,29534,24297,26078,24347,30270,24363,30014,24365,27966,24373,27454,24405,27326,24467,29310,24473,26238,24517,26878,24529,25342,24545,25086,25079,31713,25215,32657,25311,32465,25327,32209,25335,31697,25531,30577,25551,31985,25579,30193,25581,28145,25727,32649,25951,32425,25975,31657,26015,32361,26045,28521,26099,29673,26175,32537,26235,30617,26237,28569,26319,31961,26343,31193,26347,30169,26357,27609,26447,31929,26461,28345,26475,30137,26489,26553,26511,31865,26519,31353,26567,30969,26595,29177,26847,32453,26871,31685,27055,32101,27067,30565,27069,28517,27123,29669,27231,32405,27355,30421,27379,29653,27483,30389,27499,30133,27823,32077,27855,31949,27883,30157,27935,32301,28047,31853,28061,28269,28071,31085,28131,29165,28215,31517,28239,31901,28247,31389,28267,30109,28331,30045,28339,29533,28439,31293,28455,31037,28555,29821,28863,32579,28919,31683,29023,32419,29115,30563,29135,31971,29143,31459,29147,30435,29471,32307,29547,30131,29639,30963,29791,32395,29871,32075,29915,30411,29927,31179,29931,30155,30299,30363,30311,31131,30479,31803,30599,30843,30815,32391,30831,32135,31031,31527,31055,31911,31127,31335,31503,31799,31775,32271],
        },
    16 :
        {
            4 : [32790,46080,32796,39936,32809,51712,32976,34176,33034,43072,33036,38976,33042,42048,33090,41280,33096,35136,33104,34112,33300,37920,33345,49440,33348,37152,33352,35104,33376,33568,33798,45072,33840,34320,34821,53256,34849,49672,34881,49416,34882,41224,35332,36904,36994,41092,41217,49218,43009,49162,49161,51201],
            6 : [32799,64512,32862,48384,32875,60160,32919,62592,32926,48256,32935,62080,32942,47744,32971,59776,32982,46464,32995,58240,33047,62528,33070,47680,33078,46656,33138,42816,33166,47296,33189,53952,33204,38592,33208,36544,33219,57792,33222,45504,33233,50624,33323,59936,33331,58912,33433,52384,33443,58016,33460,38560,33475,57760,33505,50080,33573,53856,33577,51808,33605,53600,33633,50016,33729,49632,33898,43792,33905,50960,33912,36624,33948,40080,33989,53648,34002,42384,34017,50064,34018,41872,34146,41808,34185,51408,34186,43216,34216,35536,34311,61488,34318,47152,34356,38448,34360,36400,34374,45360,34378,43312,34385,50480,34465,49840,34497,49584,34594,41584,34626,41328,34690,41200,34692,37104,34855,61960,34862,47624,34891,59656,34958,47240,35044,37768,35113,51784,35121,50760,35141,53576,35170,41800,35341,55336,35366,45608,35395,57640,35401,51496,35428,37672,35432,35624,35459,57512,35617,49768,35910,45336,35913,51480,35924,38168,35985,50328,35986,42136,36099,57432,36164,37208,36228,37080,36372,37944,36418,41272,36894,48132,36935,61700,36941,55556,36956,40196,36963,58116,37059,57732,37066,43396,37076,38276,37139,58436,37146,44100,37161,51780,37189,53572,37193,51524,37268,38084,37282,41668,37383,61476,37398,46116,37404,39972,37413,53796,37452,39204,37514,43172,37635,57444,37642,43108,37649,50276,37697,49508,37902,47124,37972,38164,37986,41748,38021,53396,38028,39060,38049,49812,38050,41620,38081,49556,38162,42068,38210,41300,38923,59404,38931,58380,39009,49932,39043,57484,39046,45196,39298,41164,39434,43052,39442,42028,39554,41132,39942,45084,41018,44546,41038,47362,41062,45826,41125,53890,41129,51842,41161,51586,41253,53826,41349,53442,41377,49858,41409,49602,41479,61474,41507,57890,41546,43298,41553,50466,41570,41762,41605,53410,41998,47122,42010,44050,42081,49938,42117,53394,42122,43154,42249,51282,42250,43090,43015,61450,43139,57482,43269,53322,43585,49450,44035,57370,45105,50694,45189,53382,45217,49798,45321,51270,45571,57382,45601,49702,46085,53270,46593,49206,49291,59521,49301,54401,49421,55361,49477,53569,49481,51521,49539,57537,49677,55329,49689,52257,49801,51361,49925,53345,50201,52241,50225,50705,51207,61449,51213,55305,51219,58377,51235,57865,52229,53273,53379,57477],
            8 : [32991,64896,33231,63936,33262,48064,33276,40896,33375,64800,33518,48032,33525,55200,33717,55008,33718,46816,33724,40672,33784,36832,33887,64784,33903,64272,33917,57104,33982,48784,34007,62864,34035,59280,34041,53136,34110,48720,34141,56656,34222,47824,34278,46032,34281,52176,34290,42960,34365,56880,34391,62768,34425,53040,34503,61872,34537,52144,34546,42928,34638,47472,34643,58736,34701,55536,34740,38640,34755,57840,34762,43504,34776,36336,34879,65032,34999,63112,35053,56200,35065,53128,35131,61000,35159,62792,35181,56136,35190,46920,35193,53064,35229,56520,35230,48328,35258,44744,35260,40648,35278,47560,35415,62760,35446,46888,35452,40744,35527,61864,35541,54696,35562,43944,35599,63592,35611,60520,35614,48232,35661,55656,35669,54632,35692,39784,35747,58088,35753,51944,35788,39400,35793,50664,35796,38376,35895,63000,35933,56600,35964,40728,36021,54936,36026,44696,36028,40600,36057,52632,36082,42904,36139,59992,36149,54872,36209,51032,36210,42840,36245,54488,36259,58072,36261,53976,36321,50136,36430,47416,36442,44344,36453,54072,36454,45880,36468,38712,36472,36664,36491,59576,36521,51896,36524,39608,36530,42680,36547,57784,36684,39288,36705,50040,36756,38136,36801,49656,37207,62788,37291,60100,37302,46788,37361,51140,37455,63780,37491,59172,37547,60068,37550,47780,37561,52900,37564,40612,37581,55716,37618,42916,37662,48228,37788,40164,37804,39652,37809,50916,37827,57828,37842,42468,37858,41956,37949,56852,38010,44820,38094,47508,38108,40340,38115,58260,38122,43924,38173,56404,38227,58708,38246,45908,38252,39764,38291,58580,38313,51924,38348,39380,38415,63540,38429,56372,38439,62004,38446,47668,38453,54836,38457,52788,38460,40500,38492,40244,38513,50996,38554,44212,38578,42676,38663,61556,38691,57972,38693,53876,38694,45684,38754,41844,38817,49908,38971,60940,38999,62732,39022,47884,39030,46860,39034,44812,39097,52876,39123,58764,39125,54668,39142,45964,39148,39820,39191,62540,39195,60492,39198,48204,39211,59980,39213,55884,39226,44620,39246,47436,39258,44364,39276,39756,39309,55500,39324,40140,39394,41932,39501,55596,39566,47276,39589,53932,39633,50604,39693,55404,39706,44140,39777,50028,39811,57580,39813,53484,39841,49900,39951,63516,39963,60444,39982,47644,39990,46620,40035,58140,40037,54044,40102,45724,40105,51868,40114,42652,40205,55388,40234,43612,40242,42588,40354,41692,40483,57916,40485,53820,40521,51516,40579,57532,40585,51388,40593,50364,40609,49852,40641,49596,40737,49788,40738,41596,41119,64642,41175,62850,41181,56706,41182,48514,41203,59266,41209,53122,41319,62274,41334,46914,41531,60962,41551,63778,41565,56610,41627,60578,41657,52898,41671,61858,41677,55714,41678,47522,41690,44450,41774,47714,41785,52834,41829,54114,41834,43874,41929,51682,41953,50146,42015,64530,42093,56082,42135,62610,42226,42898,42263,62546,42294,46674,42298,44626,42326,46418,42330,44370,42341,54098,42406,45778,42417,50898,42511,63538,42547,58930,42553,52786,42579,58674,42595,58162,42631,61618,42662,45746,42693,53682,42694,45490,42789,53874,42819,57714,42833,50546,43099,60682,43101,56586,43115,60170,43125,55050,43182,47754,43190,46730,43219,58762,43241,52106,43310,47690,43339,59722,43411,58570,43441,50890,43563,59946,43565,55850,43610,44330,43667,58538,43670,46250,43685,53930,43721,51626,43797,54378,43798,46186,43811,57962,44062,48154,44121,52506,44137,51994,44181,54426,44197,53914,44198,45722,44230,45466,44233,51610,44241,50586,44625,50490,44675,57530,44805,53370,45115,60934,45135,63750,45166,47878,45214,48262,45235,59014,45270,46470,45289,52102,45341,56390,45365,54854,45369,52806,45383,61766,45387,59718,45462,46278,45465,52422,45481,51910,45583,63526,45597,56358,45643,59686,45645,55590,45703,61606,45717,54438,45769,51622,45955,57574,46095,63510,46158,47382,46181,54038,46185,51990,46219,59542,46243,58006,46481,50390,46606,47158,46665,51510,46753,49846,47147,59918,47149,55822,47157,54798,47182,47374,47193,52494,47329,50062,47409,50766,47441,50510,47683,57646,47761,50350,47809,49582,47905,49774,49275,61185,49371,60801,49487,63809,49501,56641,49529,53057,49579,60097,49719,63009,49751,62753,49767,62241,49819,60577,49891,58273,49947,60513,50005,54625,50055,61665,50073,52449,50089,51937,50223,64017,50327,62609,50381,55697,50387,58769,50483,58961,50489,52817,50609,50897,50763,59697,50789,54065,50835,58545,50979,57969,51315,59145,51343,63625,51381,54921,51405,55689,51539,58697,51595,59593,51603,58569,51741,56361,51769,52777,51795,58665,51853,55465,51861,54441,51875,58025,51909,53673,52247,62489,52359,61593,52373,54425,52387,58009,52421,53657,52487,61529,52995,57465,53303,62981,53309,56837,53351,62213,53403,60549,53533,56389,53645,55493,53653,54469,53799,61989,53811,58917,53861,54053,53899,59557,53907,58533,53909,54437,54301,56341,54803,58421,55325,56333,55341,55821,55373,55565,55395,58125,55619,57677,55815,61485,56835,57405,57447,62211,57499,60547,57651,58947,58395,60435,58419,58899,59527,61579,60423,61467],
            10 : [34747,61168,34791,62448,34797,56304,35581,57256,35710,49000,35743,64744,35819,60392,35967,65304,36094,49048,36311,62936,36317,56792,36318,48600,36331,60376,36471,63288,36511,64696,36527,64184,36539,61112,36639,64632,36655,64120,36669,56952,36687,63864,36711,62328,36723,59256,36729,53112,36765,56568,36779,60152,36794,44792,36813,55800,36828,40440,37503,65316,37695,65124,37757,57188,37878,47076,38207,65108,38333,57044,38334,48852,38381,56276,38523,61236,38525,57140,38637,56244,38780,40820,38814,48372,38829,56052,38835,59124,38867,58868,38889,52212,39287,63308,39375,63948,39417,53196,39550,48940,39614,48812,39659,60332,39735,63084,39798,46956,39910,46060,39921,51180,40047,64284,40061,57116,40095,64668,40111,64156,40155,60828,40158,48540,40251,61020,40286,48476,40295,62300,40313,53084,40365,56028,40398,47580,40421,54236,40425,52188,40503,63036,40606,48316,40627,59068,40629,54972,40665,52668,40678,46012,40747,60028,40775,61820,40809,52092,40846,47356,40914,42492,40930,41980,41215,65410,41407,65218,41439,64962,41463,63426,41469,57282,41723,61346,41726,49058,41839,64354,41854,48994,41918,48866,42111,65298,42223,64402,42237,57234,42351,64338,42363,61266,42365,57170,42430,48850,42459,60882,42471,62418,42478,48082,42489,53202,42619,61234,42622,48946,42671,64178,42739,59314,42741,55218,42745,53170,42807,63090,42839,62834,42869,55154,42874,44914,42923,60146,42925,56050,42934,46834,42951,61938,42993,51186,43247,64394,43451,61130,43583,65066,43631,64298,43735,62890,43770,44970,43855,63850,43898,44906,43989,54762,44143,64282,44239,63898,44263,62362,44267,60314,44350,48730,44379,60762,44395,60250,44455,62170,44459,60122,44501,54746,44515,58330,44661,55098,44699,60602,44701,56506,44757,54714,44771,58298,44827,60538,44843,60026,44851,59002,44885,54650,44998,45562,45247,65158,45501,57030,45558,47046,45691,61222,45727,64678,45811,59302,45817,53158,45917,56678,45931,60262,46005,55014,46295,62870,46311,62358,46315,60310,46395,61014,46451,59222,46507,60118,46517,54998,46535,61910,46542,47574,46549,54742,46563,58326,46647,63030,46707,59190,46749,56502,46798,47542,46809,52662,46894,47734,46899,58998,46923,59766,46925,55670,47011,58102,47230,48910,47279,64142,47287,63118,47342,48014,47407,64078,47447,62798,47463,62286,47470,47950,47565,55758,47678,48686,47707,60718,47723,60206,47821,55726,47827,58798,47901,56430,47902,48238,47915,60014,47971,58222,48013,55534,48037,53998,48159,64542,48237,56094,48245,55070,48331,59806,48399,63582,48467,58718,48469,54622,48561,50910,48667,60478,48679,62014,48741,54078,48745,52030,48803,58046,48835,57790,48917,54398,48933,53886,48965,53630,48969,51582,49029,53502,49535,65345,49631,64961,49647,64449,49855,65185,49917,57249,50043,61281,50135,62945,50367,65169,50423,63377,50429,57233,50651,60881,50669,56273,50807,63281,50931,59313,51115,60145,51177,52209,51663,63945,51693,56265,51701,55241,51705,53193,52015,64105,52027,61033,52059,60777,52077,56169,52111,63721,52181,54761,52535,63065,52587,60249,52635,60633,52653,56025,52709,54233,52935,61881,53031,62073,53035,60025,53075,58745,53091,58233,53599,64837,53735,62405,53855,64805,53979,60837,53981,56741,54107,60773,54159,63717,54195,59109,54219,59877,54383,64277,54487,62869,54515,59285,54587,61013,54607,63829,54671,63701,54699,60117,54739,58837,54831,64053,54887,62261,54951,62133,54957,55989,55115,59765,55117,55669,55359,65037,55415,63245,55483,61069,55517,56717,55539,59277,55661,56141,55719,62157,55751,61901,55863,63021,55965,56493,56035,58285,56103,62061,56435,59165,56471,62621,56619,59997,56651,59741,56723,58589,56883,58941,57099,59517,57711,64323,57807,63939,57819,60867,58043,61091,58075,60835,58479,64275,58543,64147,58551,63123,58739,59219,58767,63699,58935,63027,59035,60595,59051,60083,59179,60019,59599,63883,59627,60299,60055,62635,60559,63643,60583,62107,60807,61659,60951,62523,61503,65031,61551,64263,61743,64071,62231,62567,62519,62999,62735,63575,63535,64015],
        }
    }

Ns = [1024, 2048, 4096, 8192]

func = {
    30 : combine32_30,
    28 : combine32_28,
    26 : combine32_26,
    24 : combine32_24,
    22 : combine32_22,
    20 : combine32_20,
    18 : combine32_18
}

cursor = db.cursor()

sql = "INSERT INTO `msalgo`.`mstest_copy1` (`mstype`, `seccess`, `L`, `t`, `N`, `p`, `exectime`) VALUES ('%s', %s, %s, %s, %s, %s, %s)"

# def combine(x1,x2,x3):
#     return (x1*x2)^(x2*x3)^(x1*x3)
# z=[]
# L = 48
# p=0.75
# seed = [random.randint(0, (1 << L) - 1) for i in range(3)]
# l1 = LFSR(seed[0],0b100000000000000000000000010000000000000000000000,48)
# l2 = LFSR(seed[1],0b100000000000000000000000000000000010000000000000,48)
# l3 = LFSR(seed[2],0b100000100000000000000000000000000000000000000000,48)
# for i in range(32768):
#     z.append(combine(l1.next(), l2.next(), l3.next()))
# for N in Ns:
#     ms = MS(0b100000000000000000000000010000000000000000000000, 48, z[:N], p)
#     #print(bin(seed[0]))
#     timeA = time.perf_counter()
#     ansA = ms.crackA()
#     timeA = time.perf_counter() - timeA
#     #print(bin(ansA))
#     timeB = time.perf_counter()
#     ansB = ms.crackB()
#     timeB = time.perf_counter() - timeB
#     #print(bin(ansB))
#     if ansA == seed[0]:
#         print("A seccess")
#         cursor.execute(sql % ('A', 1, L, 2, N, p, timeA))
#     else:
#         print("A fail")
#         cursor.execute(sql % ('A', 0, L, 2, N, p, timeA))
#     if ansB == seed[0]:
#         print("B seccess")
#         cursor.execute(sql % ('B', 1, L, 2, N, p, timeB))
#     else:
#         print("B fail")
#         cursor.execute(sql % ('B', 0, L, 2, N, p, timeB))

# @timeout_decorator.timeout(600)
@func_timeout.func_set_timeout(600)
def testA(ms):
    timeA = time.perf_counter()
    z = ms.crackA()
    timeA = time.perf_counter() - timeA
    return z, timeA

# @timeout_decorator.timeout(600)
@func_timeout.func_set_timeout(600)
def testB(ms):
    timeB = time.perf_counter()
    z = ms.crackB()
    timeB = time.perf_counter() - timeB
    return z, timeB

for roundTimes in range(3):
    for correlation, combine in func.items():
        correlation = correlation / 32
        for L, diction in data.items():
            for t, maskList in diction.items():
                masks = random.sample(maskList, 5)
                seed = [random.randint(0, (1 << L) - 1) for i in range(5)]
                l1, l2, l3, l4, l5 = [LFSR(seed[i], masks[i], L) for i in range(5)]
                z = []
                for i in range(32768):
                    z.append(combine(l1.next(), l2.next(), l3.next(), l4.next(), l5.next()))
                for N in Ns:
                    ms = MS(masks[0], L, z[:N], correlation)
                    try:
                        ansA, timeA = testA(ms)
                    except func_timeout.exceptions.FunctionTimedOut:
                        print("A timeout")
                        cursor.execute(sql % ('A', 0, L, t, N, correlation, 'null'))
                    else:
                        if ansA == seed[0]:
                            print("A seccess")
                            cursor.execute(sql % ('A', 1, L, t, N, correlation, timeA))
                        else:
                            print("A fail")
                            cursor.execute(sql % ('A', 0, L, t, N, correlation, timeA))
                    try:
                        ansB, timeB = testB(ms)
                    except func_timeout.exceptions.FunctionTimedOut:
                        print("B timeout")
                        cursor.execute(sql % ('B', 0, L, t, N, correlation, 'null'))
                    else:
                        if ansB == seed[0]:
                            print("B seccess")
                            cursor.execute(sql % ('B', 1, L, t, N, correlation, timeB))
                        else:
                            print("B fail")
                            cursor.execute(sql % ('B', 0, L, t, N, correlation, timeB))

                    # timeA = time.perf_counter()
                    # ansA = ms.crackA()
                    # timeA = time.perf_counter() - timeA
                    # if ansA == seed[0]:
                    #     print("A seccess")
                    #     cursor.execute(sql % ('A', 1, L, t, N, correlation, timeA))
                    # else:
                    #     print("A fail")
                    #     cursor.execute(sql % ('A', 0, L, t, N, correlation, timeA))
                    # timeB = time.perf_counter()
                    # ansB = ms.crackB()
                    # timeB = time.perf_counter() - timeB
                    # if ansB == seed[0]:
                    #     print("B seccess")
                    #     cursor.execute(sql % ('B', 1, L, t, N, correlation, timeB))
                    # else:
                    #     print("B fail")
                    #     cursor.execute(sql % ('B', 0, L, t, N, correlation, timeB))
                    db.commit()

db.close()