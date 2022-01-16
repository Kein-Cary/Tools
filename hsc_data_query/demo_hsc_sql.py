sql_str = """
    SELECT
        object_id, ra, dec,
        a_g, a_r, a_i,
        g_cmodel_mag, r_cmodel_mag, i_cmodel_mag,
        specz_redshift

    FROM
        pdr3_wide.forced
        LEFT JOIN pdr3_wide.specz AS specz USING (object_id)
    WHERE
        isprimary
        AND g_cmodel_mag >0
            AND NOT g_cmodel_mag = 'NaN'
            AND NOT g_cmodel_mag = 'inf'
        AND r_cmodel_mag >0
            AND NOT r_cmodel_mag = 'NaN'
            AND NOT r_cmodel_mag = 'inf'
        AND i_cmodel_mag >0
            AND NOT i_cmodel_mag = 'NaN'
            AND NOT i_cmodel_mag = 'inf'

        AND i_extendedness_value = 1

        AND NOT i_extendedness_flag
        AND NOT g_cmodel_flag
        AND NOT r_cmodel_flag
        AND NOT i_cmodel_flag
        AND i_cmodel_mag >= mag_0
        AND i_cmodel_mag <= mag_1
"""

sql_str = """
    WITH
        -- Describe the catalog to match
        my_catalog( my_id, my_ra, my_dec ) AS (VALUES
            (NULL::int8, NULL::float8, NULL::float8 ),
            (43224722665657008, 172.46470018, 2.93796983),
            (44191588523470528, 140.53187849, 3.76633851),
            (43251360052824696, 180.97275624, 1.78346317),
            (43211240763322376, 168.19808871, 2.50297013),
            (41693914716982864, 16.54324276, 0.85571472),
            (44240667114751600, 155.91947755, 4.18753349),
            (44777344753224256, 337.36449932, 3.28870717),
            (71178671930498040, 245.22622775, 44.70047377),
            (41663257240422056, 6.00128458, 0.50204562),
            (70373155109097920, 228.79440527, 43.66433032),
            (43752329333208936, 351.09058477, 2.10517341),
            (40124078400483032, 206.56875652, -2.46142333),
            (42287376118087976, 216.84119854, 0.94553458),
            (41637023580188784, 356.86499314, -0.15380983),
            (43237607567550856, 177.26478023, 1.50174496),
            (70391434489918416, 235.58189731, 43.75787164),
            (41646094551111960, 359.46881319, 1.27529443),
            (44271307411450616, 166.52467523, 3.98917648),
            (45823963858748784, 331.52825896, 5.68886881),
            (44201355279095584, 142.25656562, 4.09368945),
            (44429770229830080, 220.06090113, 3.69501243),
            (42269929960921024, 210.69993979, 1.22207444),
            (45261404747356048, 140.63820955, 4.81712532),
            (41734046891393424, 29.16192644, 0.84213809),
            (42727034740297136, 5.20113840, 2.04945595),
            (43348254515026032, 213.52599572, 1.87023132),
            (42247797994446072, 203.50708774, 1.19982188),
            (70351989510266120, 217.15937177, 43.64412001),
            (41693760098161600, 16.62842146, 0.13723751),
            (44248913451961288, 159.66974836, 4.23135848),
            (45256319506075056, 140.09080176, 4.83973123),
            (40677948793047400, 34.20439497, -1.34814834),
            (70409155524978200, 243.44396725, 43.26837972),
            (39622705393197336, 37.01131655, -2.20050165),
            (42256710051584240, 206.34469014, 0.20602360),
            (42797691247280136, 28.72023958, 2.68870557),
            (69634708497004064, 248.78587655, 42.80004536),
            (43286282431905416, 193.24545805, 2.32667317),
            (39614308732131216, 33.50057683, -2.84845083),
            (69581923348927192, 224.87788081, 42.38190604),
            (45854866148446240, 341.76590593, 4.77395408),
            (40594411679141632, 6.14677609, -0.21866742),
            (45855020767271216, 341.65737162, 5.65967961),
            (41671915894490144, 9.06757449, 0.46398513),
            (45827674710493560, 333.92394801, 5.68815490),
            (43312657826070880, 202.19793836, 1.76594970),
            (42288351075654384, 215.62197182, 1.43498094),
            (41228005254650936, 219.05686947, -0.26934002),
            (38562767594081152, 38.41396969, -4.17528446),
            (40625734875628328, 15.65211760, -0.83547062),
            (45850210403897376, 340.58547510, 5.55958760),
            (43250810297012304, 181.62657627, 1.82071768),
            (40625468587659400, 16.10340325, -0.41418467),
            (40598792545782056, 7.63105200, -0.93126505),
            (69581931938861752, 224.93369027, 42.71229832),
            (44781077079806872, 339.69738966, 4.24750001),
            (43356217384392560, 217.57930569, 1.49321929),
            (45283257540958432, 148.18571632, 4.88191512),
            (38540923390418168, 30.85046815, -3.76801925),
            (44438699466841008, 222.77263026, 3.53702290),
            (45260588703562896, 141.79058833, 5.21054220),
            (43264700221253480, 185.21493730, 2.27029366),
            (41671907304559312, 9.06174079, 0.16114470),
            (43303582560181296, 199.55917739, 1.71328507),
            (69612585120457184, 239.12895763, 42.87910659),
            (37476716393813936, 30.59029984, -5.90236492),
            (44196948642648464, 140.68649224, 3.70678801),
            (41165865667810672, 199.05339377, -0.93369787),
            (41592497654223136, 342.84785609, -0.00695634),
            (70413012405610392, 246.50846387, 43.70414631),
            (43225122097619672, 171.89483845, 2.47857387),
            (43765815530506640, 355.24776190, 2.74564764),
            (43264000141580720, 186.10274194, 1.73334841),
            (42771994457947120, 18.86304246, 2.79782738),
            (43149947285040600, 147.00473673, 2.68082543),
            (41645531910399496, 0.21053670, 0.78477247),
            (40564316843291120, 354.84035734, -1.50707774),
            (41667801315823240, 7.30240665, 0.97226510),
            (46892556016954656, 334.51780884, 7.41042783),
            (42758100238751368, 15.33873030, 2.29891400),
            (44222791460875616, 150.35981724, 3.98505835),
            (43145686677479488, 145.42146046, 2.66799687),
            (41671769865607136, 9.20060983, 0.15762638),
            (69621793530337664, 242.33905293, 42.91910813),
            (44293709960863104, 173.45089514, 3.94285400),
            (70360111293421432, 222.43406874, 44.18024731),
            (43347979637117616, 213.86403044, 1.85057152),
            (42644167641289256, 336.05676866, 0.84155046),
            (40171786897215072, 223.81311915, -1.62911962),
            (42780760486205440, 21.81707106, 1.54420603),
            (40625593141708888, 15.96992975, -0.96044457),
            (43255779574177896, 182.33701455, 2.82213705),
            (42705598558531776, 357.06320732, 0.78526342),
            (38545454580910464, 32.18810385, -3.97824274),
            (41746966153023344, 33.96368143, 0.87012302),
            (69586729417337848, 226.26975956, 42.27591933),
            (39622975976132912, 36.58471166, -2.47337862),
            (42635225519384912, 333.36695233, 0.54423680),
            (43717157846018176, 339.12543297, 2.60838823),
            (42287517852007392, 216.63163613, 1.08757686),
            (37471644037433728, 30.13095118, -5.46698663),
            (46888136495615456, 333.11677028, 6.49194447),
            (41632904706543392, 355.07289622, -0.06583638),
            (39613501278281920, 34.59952780, -1.98425819),
            (44240113063976336, 156.72944637, 4.14563558),
            (42300823660685152, 220.86125879, 0.01129965),
            (39600040850773976, 30.44209285, -1.64288570),
            (42643768209333584, 336.69115121, 1.41823316),
            (39604718070155664, 31.56600147, -1.51276394),
            (45317900747162896, 160.84095821, 5.21196999),
            (44315433905443336, 181.21334353, 4.21650085),
            (42688298430261736, 350.71602361, 1.46329793),
            (41610201509414816, 348.57408030, -1.09511868),
            (41218796844784480, 216.63912155, -0.10900476),
            (69559658238464064, 215.48135759, 42.34048797),
            (44235556103669768, 155.44025623, 3.12991978),
            (44438579207744576, 223.03224786, 4.21697978),
            (42287371823120720, 216.77750064, 0.72077356),
            (41747378469885296, 33.36928694, 1.00076450),
            (40154740172019000, 217.09736099, -1.80832152),
            (38562235018143256, 39.26355300, -3.29336193),
            (42279121190943448, 213.09927188, 0.63970972),
            (40681938817671616, 36.29468729, -1.08316488),
            (42763172595125728, 15.87139559, 1.66336089),
            (43145828411392352, 145.25875960, 2.69901738),
            (42283261539413688, 214.96281736, 1.32009940),
            (70417260128267024, 248.86293886, 43.16392019),
            (42301403481268920, 220.22003533, 1.31362410),
            (39627112029638216, 38.48387824, -1.88464676),
            (44240366467042040, 156.32969665, 3.12353150),
            (42234998991905536, 198.47058278, 0.40438482),
            (69569145821222592, 218.26868952, 42.55033802),
            (41676459969893296, 10.31091525, 0.98544322),
            (42700917044170912, 355.89860154, 0.33094159),
            (42643897058355552, 336.42815039, 1.06081107),
            (69634991964843616, 248.48235095, 43.09116930),
            (42269929960927096, 210.65882190, 1.35967542),
            (42309344875802936, 224.30630800, 0.08445105),
            (42829010148812960, 38.38232266, 1.91081205),
            (40594407384170688, 6.04858969, -0.44744441),
            (42833700253090176, 39.48848159, 2.57233336),
            (41759593356876280, 39.10768090, 0.23517348),
            (44429207589116952, 220.81555796, 3.20005309),
            (37498577777352600, 38.21534905, -5.56636117),
            (41694030681102544, 16.23069614, 0.06009439),
            (43264283609424376, 185.67209979, 2.09973340),
            (41729236528028224, 28.17483344, 1.00708233),
            (43774040392879608, 358.91159770, 1.76528705),
            (40677695389977216, 34.62748105, -0.42965444),
            (44209721875387088, 145.71205205, 3.32135857),
            (42705873436436152, 356.59953702, 0.74942440),
            (37480435835489456, 33.11248520, -5.62656508),
            (69603647293518040, 235.48775148, 42.82190598),
            (43145141216625640, 146.19638271, 2.73271243),
            (44442680901525536, 224.86562265, 3.34607977),
            (44790702101512960, 341.60027980, 4.38386823),
            (41601272272415504, 345.81976221, -0.85392949),
            (69586467424329216, 226.84656313, 42.73408072),
            (36420244633308336, 33.60742860, -6.46116149),
            (70412741822667048, 247.06315244, 43.81329997),
            (44240525380831568, 156.25668630, 4.03402979),
            (43774182126801984, 358.78608003, 1.97939668),
            (42300690516699488, 221.03777056, 0.17847792),
            (42798077794339400, 28.08407179, 1.55732689),
            (44350884565506736, 192.87551423, 3.96769598),
            (41672049038481616, 8.97972971, 0.35647855),
            (70373292548049408, 228.59714125, 43.64938659),
            (70365041915877344, 223.36648016, 43.49660821),
            (42287792729913912, 216.24291275, 1.15360899),
            (42710542065887696, 357.83886952, 0.61693608),
            (43760584260343416, 354.75337408, 2.51097448),
            (43268969418744816, 186.82361894, 2.60164577),
            (37476170932967768, 31.36762409, -5.73242647),
            (45278322623527360, 147.55058285, 5.43814046),
            (43765952969458400, 354.94375996, 2.77538625),
            (42758087353847888, 15.23449645, 1.66606345),
            (42718517820150992, 1.83321195, 2.25365486),
            (44201497013022448, 142.08658523, 4.36230359),
            (41694314148940640, 15.84749801, 0.31946317),
            (44231849546902952, 153.08394355, 3.42742784),
            (42710404626932496, 357.93294402, 0.56080350),
            (37485383637817048, 33.87296546, -5.54867276),
            (41751493048550808, 35.13743693, 0.46762552),
            (40673705365356864, 32.66621403, -0.64211011),
            (45287101536684288, 150.36981474, 4.57701750),
            (69577667036339552, 222.87296613, 42.59295803),
            (70364363311045960, 224.71899448, 43.88034444),
            (43224447787751264, 172.72183042, 2.91563781),
            (41236522174803656, 222.49224019, -0.33080741),
            (45256328096004768, 140.10423746, 5.11945443),
            (41637298458089120, 356.51952476, -0.18573951),
            (43215217903037696, 170.24855267, 2.13785592),
            (40599496920415968, 6.70326004, -0.26935616),
            (45260438379716512, 141.86822776, 4.66959719),
            (42247944023326976, 203.33666933, 1.41826907),
            (42793958920696944, 26.30360819, 1.61953892),
            (70356524995728456, 218.93117600, 43.62635926),
            (43219895122417528, 171.50212237, 2.16219171),
            (40630283246003040, 17.08694291, -0.15459299),
            (42243382768066816, 202.08403337, 0.43852623),
            (40599226337478832, 7.08764651, -0.01709668),
            (44275851486844192, 167.90486300, 4.21968034),
            (45260584408597448, 141.76883162, 4.98207173),
            (37498732396177776, 38.13852219, -4.77523519),
            (44254140427160168, 160.14928952, 4.40868061),
            (37484829587038416, 34.65519704, -5.67358956),
            (41637276983246928, 356.45069783, -1.12774958),
            (42239680506250512, 199.69731178, 0.68020649),
            (41671628131689680, 9.41911426, 0.08800863),
            (39600848304627432, 29.33376998, -2.40237388),
            (43224851514676496, 172.27894833, 2.55202692),
            (42775696719765288, 21.22434214, 2.50825399),
            (41170693211047496, 200.04549960, -0.26899014),
            (41750823033660800, 36.18028800, 1.38007448),
            (70391705072853232, 235.27065591, 43.45078012),
            (41693622659208080, 16.86172222, 0.14714311),
            (42257002109362696, 205.94324351, 0.97912114),
            (44213857928895728, 147.73094918, 3.98923281),
            (41165861372844496, 199.13054889, -1.08661843),
            (44231596143829880, 153.48032252, 4.33128034),
            (69613392574307000, 237.73621427, 42.24406374),
            (42287934463829848, 216.19778561, 1.33630104),
            (40633706334932456, 19.85983353, -0.74364964),
            (45278026270797896, 147.83728121, 4.54262518),
            (40101413858068872, 200.07996964, -1.83245064),
            (38545055148961808, 32.75180923, -3.28503947),
            (69625766375086576, 245.15506273, 42.38081438),
            (69603926466390088, 235.07211956, 43.00676931),
            (45318038186115584, 160.71644149, 5.18920125),
            (41575021432298144, 336.73221565, -0.91765984),
            (44350613982561776, 193.20086657, 4.02235416),
            (42692262685077448, 352.86019312, 0.61524919),
            (40598521962847752, 8.04672486, -0.66670415),
            (44333292379460064, 186.79872461, 3.86699531),
            (44776945321267616, 337.96065946, 3.90808996),
            (43329970839242600, 208.51164852, 1.60830705),
            (42648415363948400, 337.85897917, 0.29530483),
            (43260284994872448, 183.74275216, 1.51130676),
            (40608151279516912, 9.73852647, -0.44944086),
            (70373288253088560, 228.48806459, 43.55823624),
            (69625646115999936, 245.46404612, 43.09366016),
            (45265545095816032, 142.47410114, 5.53441112),
            (39609635807716240, 32.42673643, -2.75061308),
            (41227854930801648, 219.21215525, -0.70025940),
            (41597024549753032, 344.16753378, -0.41416100),
            (41712069543741520, 21.65829397, 1.42999763),
            (70347716017807776, 214.91311858, 43.19543431),
            (69595388071397872, 230.59024961, 42.20478705),
            (41746811534207352, 34.08554275, 0.26577317),
            (70382234669964496, 232.30930661, 44.03101967),
            (41228138398643824, 218.79727598, -0.29333685),
            (44205336713774496, 144.28766445, 3.81442656),
            (70387178177322288, 233.32958708, 43.85949915),
            (37489639950405192, 35.48496093, -5.75879876),
            (40625335443673920, 16.29322940, -0.15980839),
            (70355837800958600, 220.26009112, 43.61500588),
            (45854200428516992, 342.73764357, 5.82523945),
            (44262086116662728, 164.21622516, 3.29462735),
            (45304569168670528, 156.52757007, 5.14354761),
            (69585900488636472, 227.72139834, 42.10852346),
            (42644021612409232, 336.26016059, 0.55177412),
            (41734034006496216, 29.09469014, 0.41038466),
            (38545196882878544, 32.57302588, -3.16523129),
            (38549457490437352, 34.27636846, -3.09662239),
            (69630851616368248, 245.81807767, 42.43516271),
            (41165998811798552, 198.91593950, -1.07465256),
            (42727447057157200, 4.62013063, 2.01291632),
            (44210001048258864, 145.32256979, 3.46678450),
            (44205070425804128, 144.63474048, 4.25187813),
            (42731557340860336, 6.45815899, 1.54972464),
            (41698725080358616, 17.34616958, 0.90818140),
            (42683616915903328, 349.62664405, 1.09202851),
            (44777357638126768, 337.39390059, 3.84343094),
            (41188676239126592, 205.35591775, -1.04186332),
            (42779948737381560, 22.81438240, 2.08244542),
            (43260014411933896, 184.16821401, 1.70164122),
            (41676867991784248, 9.82501063, 0.69983168),
            (40110888555920768, 202.12910450, -2.19736274),
            (41720586463885552, 24.90972400, 1.31854015),
            (37497899172517728, 39.20028704, -5.20407545),
            (42692400124028336, 352.66733942, 0.60938179),
            (41698729375320496, 17.36377615, 1.03161492),
            (69581661355923560, 225.63918494, 42.93234679),
            (44270757655630496, 167.28900849, 3.82228024),
            (41658433992152992, 5.06713357, 0.07947047),
            (43127549030587360, 140.20710220, 2.75388558),
            (41663798406303264, 5.20390356, 0.18223493),
            (44390174926328392, 206.51473008, 3.15860901),
            (40665025236457800, 29.36350555, -1.39712263),
            (43145819821464528, 145.21681918, 2.45462554),
            (42701458210056816, 355.24903202, 0.08159323),
            (38567577957451240, 39.46805138, -4.18195553),
            (46892951153961336, 333.99885575, 6.73070037),
            (43145957260414576, 145.07137186, 2.41798703),
            (41219338010655440, 215.83133253, -0.53893477),
            (44205757620574320, 143.70457749, 4.24452603),
            (73979707867024400, 213.62542588, 51.93785185),
            (41566362778226552, 333.60311748, -0.91005711),
            (42766621453865360, 18.57829373, 2.25572801),
            (40665171265339248, 29.18151946, -1.19392862),
            (40120101260768560, 204.57415928, -2.10707183),
            (41746399217337528, 34.67478439, 0.11377806),
            (40150621298384056, 215.29993131, -1.64552901),
            (42662026115310208, 341.63083394, 0.58320984),
            (40629600346200520, 17.95341745, -0.01815485),
            (43127394411759568, 140.39419476, 1.97325867)
        )
        ,

        match AS (
            SELECT
                my_catalog.*,
                object_id

            FROM
                my_catalog JOIN pdr3_wide.forced ON (
                -- Match objects within 1.0 arcsec
                coneSearch( coord, my_ra, my_dec, 1.0 )
                )
            )
    -- The principal clause:
    --   Filter the above match catalog by conditions
    --   and select columns other than object_id
    SELECT
        match.*,
        object_id,
        g_apertureflux_10_flux, r_apertureflux_10_flux, i_apertureflux_10_flux,
        g_apertureflux_10_fluxerr, r_apertureflux_10_fluxerr, i_apertureflux_10_fluxerr,

        g_apertureflux_15_flux, r_apertureflux_15_flux, i_apertureflux_15_flux,
        g_apertureflux_15_fluxerr, r_apertureflux_15_fluxerr, i_apertureflux_15_fluxerr,

        g_apertureflux_20_flux, r_apertureflux_20_flux, i_apertureflux_20_flux,
        g_apertureflux_20_fluxerr, r_apertureflux_20_fluxerr, i_apertureflux_20_fluxerr,

        g_apertureflux_30_flux, r_apertureflux_30_flux, i_apertureflux_30_flux,
        g_apertureflux_30_fluxerr, r_apertureflux_30_fluxerr, i_apertureflux_30_fluxerr,

        g_apertureflux_40_flux, r_apertureflux_40_flux, i_apertureflux_40_flux,
        g_apertureflux_40_fluxerr, r_apertureflux_40_fluxerr, i_apertureflux_40_fluxerr,

        g_apertureflux_57_flux, r_apertureflux_57_flux, i_apertureflux_57_flux,
        g_apertureflux_57_fluxerr, r_apertureflux_57_fluxerr, i_apertureflux_57_fluxerr,

        g_apertureflux_84_flux, r_apertureflux_84_flux, i_apertureflux_84_flux,
        g_apertureflux_84_fluxerr, r_apertureflux_84_fluxerr, i_apertureflux_84_fluxerr,

        g_apertureflux_118_flux, r_apertureflux_118_flux, i_apertureflux_118_flux,
        g_apertureflux_118_fluxerr, r_apertureflux_118_fluxerr, i_apertureflux_118_fluxerr,

        g_apertureflux_168_flux, r_apertureflux_168_flux, i_apertureflux_168_flux,
        g_apertureflux_168_fluxerr, r_apertureflux_168_fluxerr, i_apertureflux_168_fluxerr,

        g_apertureflux_235_flux, r_apertureflux_235_flux, i_apertureflux_235_flux,
        g_apertureflux_235_fluxerr, r_apertureflux_235_fluxerr, i_apertureflux_235_fluxerr

    FROM
        match JOIN pdr3_wide.forced3 USING (object_id)
"""

"""
    -- limited catalog, so the condition is not necessary
    WHERE
        isprimary
        AND g_cmodel_mag >0
            AND NOT g_cmodel_mag = 'NaN'
            AND NOT g_cmodel_mag = 'inf'
        AND r_cmodel_mag >0
            AND NOT r_cmodel_mag = 'NaN'
            AND NOT r_cmodel_mag = 'inf'
        AND i_cmodel_mag >0
            AND NOT i_cmodel_mag = 'NaN'
            AND NOT i_cmodel_mag = 'inf'
        AND NOT g_cmodel_flag
        AND NOT r_cmodel_flag
        AND NOT i_cmodel_flag
"""