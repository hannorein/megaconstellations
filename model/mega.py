import numpy as np
import rebound
MEarth = 5.97e24
REarth = 6378.135e3

constellations_all = {
    "Sunrise": [ 
        { 'NPLANES': 1, 'ALT':500.0 ,  'INC': 97.4   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':510.3 ,  'INC': 97.4   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':520.7 ,  'INC': 97.5   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':531.0 ,  'INC': 97.5   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':541.4 ,  'INC': 97.6   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':551.7 ,  'INC': 97.6   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':562.1 ,  'INC': 97.7   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':572.4 ,  'INC': 97.7   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':582.8 ,  'INC': 97.7   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':593.1 ,  'INC': 97.8   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':603.5 ,  'INC': 97.8   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':613.8 ,  'INC': 97.9   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':624.1 ,  'INC': 97.9   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':634.5 ,  'INC': 97.9   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':644.8 ,  'INC': 98.0   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':655.2 ,  'INC': 98.0   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':665.5 ,  'INC': 98.1   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':675.9 ,  'INC': 98.1   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':686.2 ,  'INC': 98.1   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':696.5 ,  'INC': 98.2   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':706.9 ,  'INC': 98.2   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':717.2 ,  'INC': 98.3   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':727.6 ,  'INC': 98.3   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':737.9 ,  'INC': 98.3   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':748.3 ,  'INC': 98.4   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':758.6 ,  'INC': 98.4   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':769.0 ,  'INC': 98.5   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':779.3 ,  'INC': 98.5   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':789.7 ,  'INC': 98.6   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':800.0 ,  'INC': 98.6   , 'SATPP':  740 }, 
        { 'NPLANES': 1, 'ALT':810.0 ,  'INC': 98.7   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':820.2 ,  'INC': 98.7   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':830.5 ,  'INC': 98.7   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':840.7 ,  'INC': 98.8   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':850.9 ,  'INC': 98.8   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':861.2 ,  'INC': 98.9   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':871.4 ,  'INC': 98.9   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':881.6 ,  'INC': 99.0   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':891.9 ,  'INC': 99.0   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':902.1 ,  'INC': 99.1   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':912.3 ,  'INC': 99.1   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':922.6 ,  'INC': 99.2   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':932.8 ,  'INC': 99.2   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':943.0 ,  'INC': 99.3   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':953.3 ,  'INC': 99.3   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':963.5 ,  'INC': 99.3   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':973.7 ,  'INC': 99.4   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':984.0 ,  'INC': 99.4   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':994.2 ,  'INC': 99.5   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':1004.4,  'INC': 99.5   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':1014.6,  'INC': 99.6   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':1024.9,  'INC': 99.6   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':1035.1,  'INC': 99.7   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':1045.3,  'INC': 99.7   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':1055.6,  'INC': 99.8   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':1065.8,  'INC': 99.8   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':1076.0,  'INC': 99.9   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':1086.3,  'INC': 99.9   , 'SATPP':  300 }, 
        { 'NPLANES': 1, 'ALT':1096.5,  'INC': 100.0  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1106.7,  'INC': 100.0  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1117.0,  'INC': 100.0  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1127.2,  'INC': 100.1  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1137.4,  'INC': 100.2  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1147.7,  'INC': 100.2  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1157.9,  'INC': 100.2  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1168.1,  'INC': 100.3  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1178.4,  'INC': 100.3  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1188.6,  'INC': 100.4  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1198.8,  'INC': 100.4  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1209.1,  'INC': 100.5  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1219.3,  'INC': 100.5  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1229.5,  'INC': 100.6  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1239.8,  'INC': 100.6  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1250.0,  'INC': 100.7  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1260.0,  'INC': 100.7  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1270.2,  'INC': 100.8  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1280.4,  'INC': 100.8  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1290.6,  'INC': 100.9  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1300.8,  'INC': 100.9  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1310.9,  'INC': 101.0  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1321.1,  'INC': 101.0  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1331.3,  'INC': 101.1  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1341.5,  'INC': 101.2  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1351.7,  'INC': 101.2  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1361.9,  'INC': 101.3  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1372.1,  'INC': 101.3  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1382.3,  'INC': 101.4  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1392.5,  'INC': 101.4  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1402.6,  'INC': 101.5  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1412.8,  'INC': 101.5  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1423.0,  'INC': 101.6  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1433.2,  'INC': 101.7  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1443.4,  'INC': 101.7  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1453.6,  'INC': 101.8  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1463.8,  'INC': 101.8  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1474.0,  'INC': 101.9  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1484.2,  'INC': 101.9  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1494.3,  'INC': 102.0  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1504.5,  'INC': 102.0  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1514.7,  'INC': 102.1  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1524.9,  'INC': 102.2  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1535.1,  'INC': 102.2  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1545.3,  'INC': 102.3  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1555.5,  'INC': 102.3  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1565.7,  'INC': 102.4  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1575.8,  'INC': 102.4  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1586.0,  'INC': 102.5  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1596.2,  'INC': 102.5  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1606.4,  'INC': 102.6  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1616.6,  'INC': 102.7  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1626.8,  'INC': 102.7  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1637.0,  'INC': 102.8  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1647.2,  'INC': 102.8  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1657.4,  'INC': 102.9  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1667.5,  'INC': 102.9  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1677.7,  'INC': 103.0  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1687.9,  'INC': 103.0  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1698.1,  'INC': 103.1  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1708.3,  'INC': 103.2  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1718.5,  'INC': 103.2  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1728.7,  'INC': 103.3  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1738.9,  'INC': 103.3  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1749.1,  'INC': 103.4  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1759.2,  'INC': 103.4  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1769.4,  'INC': 103.5  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1779.6,  'INC': 103.5  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1789.8,  'INC': 103.6  , 'SATPP':   300}, 
        { 'NPLANES': 1, 'ALT':1800.0,  'INC': 103.7  , 'SATPP':   300}, 
              ],
    "SXODC": [
            { 'ALT':550.0  , 'INC': 32.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':552.0  , 'INC': 31.3  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':554.0  , 'INC': 30.7  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':556.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':558.0  , 'INC': 29.3  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':560.0  , 'INC': 28.7  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':562.0  , 'INC': 28.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':564.0  , 'INC': 27.3  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':566.0  , 'INC': 26.7  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':568.0  , 'INC': 26.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':565.0  , 'INC': 97.7  , 'NPLANES': 2  , 'SATPP': 4999   },
            { 'ALT':567.2  , 'INC': 97.7  , 'NPLANES': 2  , 'SATPP': 4999   },
            { 'ALT':569.4  , 'INC': 97.7  , 'NPLANES': 2  , 'SATPP': 4999   },
            { 'ALT':571.7  , 'INC': 97.7  , 'NPLANES': 2  , 'SATPP': 4999   },
            { 'ALT':573.9  , 'INC': 97.7  , 'NPLANES': 2  , 'SATPP': 4999   },
            { 'ALT':576.1  , 'INC': 97.7  , 'NPLANES': 2  , 'SATPP': 4999   },
            { 'ALT':578.3  , 'INC': 97.7  , 'NPLANES': 2  , 'SATPP': 4999   },
            { 'ALT':580.6  , 'INC': 97.7  , 'NPLANES': 2  , 'SATPP': 4999   },
            { 'ALT':582.8  , 'INC': 97.7  , 'NPLANES': 2  , 'SATPP': 4999   },
            { 'ALT':585.0  , 'INC': 97.7  , 'NPLANES': 2  , 'SATPP': 4999   },
            { 'ALT':686.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':687.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':688.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':690.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':691.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':692.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':694.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':695.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':696.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':698.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':699.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':700.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':702.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':703.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':704.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':706.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':707.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':708.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':710.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':711.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':712.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':714.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':715.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':716.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':718.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':707.0  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':708.8  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':710.5  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':712.3  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':714.0  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':715.8  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':717.6  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':719.3  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':721.1  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':722.9  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':724.6  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':726.4  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':728.1  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':729.9  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':731.7  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':733.4  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':735.2  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':737.0  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':738.7  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':740.5  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':742.2  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':744.0  , 'INC': 97.2  , 'NPLANES': 2  , 'SATPP': 5565   },
            { 'ALT':946.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':947.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':948.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':950.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':951.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':952.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':954.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':955.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':956.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':958.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':959.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':960.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':962.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':963.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':964.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':966.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':967.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':968.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':970.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':971.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':972.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':974.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':975.3  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':976.7  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':978.0  , 'INC': 30.0  , 'NPLANES': 30 , 'SATPP':  333   },
            { 'ALT':967.0  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':968.7  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':970.3  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':972.0  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':973.7  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':975.3  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':977.0  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':978.7  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':980.3  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':982.0  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':983.7  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':985.3  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':987.0  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':988.7  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':990.3  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':992.0  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':993.7  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':995.3  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':997.0  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':998.7  , 'INC': 99.4  , 'NPLANES': 2  , 'SATPP': 5770   },
            { 'ALT':1000.3 , 'INC':  99.4 , 'NPLANES':  2 , 'SATPP':  5770  },
            { 'ALT':1002.0 , 'INC':  99.4 , 'NPLANES':  2 , 'SATPP':  5770  },
              ],
    "Starlink": [ {'NPLANES':7178,'SATPP':1,'INC':30,'ALT':328},
        {'NPLANES':7178,'SATPP':1,'INC':40,'ALT':334},
        {'NPLANES':7178,'SATPP':1,'INC':53,'ALT':345},
        {'NPLANES':40,'SATPP':50,'INC':96.9,'ALT':360},
        {'NPLANES':1998,'SATPP':1,'INC':75,'ALT':373},
        {'NPLANES':4000,'SATPP':1,'INC':53,'ALT':499},
        {'NPLANES':12,'SATPP':12,'INC':148,'ALT':604},
        {'NPLANES':18,'SATPP':18,'INC':115.7,'ALT':614},
        {'NPLANES':2547,'SATPP':1,'INC':53,'ALT':345.6},
        {'NPLANES':2478,'SATPP':1,'INC':48,'ALT':340.8},
        {'NPLANES':2493,'SATPP':1,'INC':42,'ALT':335.9},
        {'NPLANES':32,'SATPP':50,'INC':53,'ALT':550},
        {'NPLANES':72,'SATPP':22,'INC':53.2,'ALT':540},
        {'NPLANES':36,'SATPP':20,'INC':70,'ALT':570},
        {'NPLANES':6,'SATPP':58,'INC':97.6,'ALT':560},
        {'NPLANES':4,'SATPP':43,'INC':97.6,'ALT':560.1},],
    "OneWeb": [ {'NPLANES':18,'SATPP':40,'INC':87.9,'ALT':1200},
        {'NPLANES':36,'SATPP':49,'INC':87.9,'ALT':1200},
        {'NPLANES':32,'SATPP':72,'INC':40,'ALT':1200},
        {'NPLANES':32,'SATPP':72,'INC':55,'ALT':1200},],
    "StarNet/GW": [ {'NPLANES':16,'SATPP':30,'INC':85,'ALT':590},
        {'NPLANES':40,'SATPP':50,'INC':50,'ALT':600},
        {'NPLANES':60,'SATPP':60,'INC':55,'ALT':508},
        {'NPLANES':48,'SATPP':36,'INC':30,'ALT':1145},
        {'NPLANES':48,'SATPP':36,'INC':40,'ALT':1145},
        {'NPLANES':48,'SATPP':36,'INC':50,'ALT':1145},
        {'NPLANES':48,'SATPP':36,'INC':60,'ALT':1145},],
    "Kuiper": [ {'NPLANES':34,'SATPP':34,'INC':51.9,'ALT':630},
        {'NPLANES':36,'SATPP':36,'INC':42,'ALT':610},
        {'NPLANES':28,'SATPP':28,'INC':33,'ALT':509},],
    }

def getAirmass(z):
    # z is the Zenith enagle
    X = 1./(np.cos(z) + 0.50572*(6.07995+90-z*180/np.pi)**(-1.6364))  # Kasten and Young (1989)
    #X = 1./np.cos(z) * (1-0.0012*np.tan(z)**2)  # Young and Irvine (1967)
    return X

def add_to_simulation(sim, ICs, debug=False):
    for IC in ICs:
        nplanes=IC['NPLANES']
        nsat=IC['SATPP']
        a = IC['ALT']*1000.+REarth
        Omegas = np.linspace(0.,2.*np.pi,nplanes)
        for i, Omega in enumerate(Omegas):
            # 5 percent jitter
            Ms = np.linspace(0.,2.*np.pi,nsat)+ 2.*np.pi/nsat*0.25*np.random.normal(size=nsat)
            for j, M in enumerate(Ms):
                sim.add(primary=sim.particles[0], M=M, a=a, omega=0, e=0, Omega=Omega, inc=IC['INC']*np.pi/180.)
                if debug and sim.N>100:
                    return
def rotY(xyz,alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    M = np.array([[c,0,-s],[0,1,0],[s,0,c]])
    return xyz @ M
def rotZ(xyz,alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    M = np.array([[c,-s,0],[s,c,0],[0,0,1]])
    return xyz @ M

def length_of_night(month,latitude, p=0):
    # https://www.ikhebeenvraag.be/mediastorage/FSDocument/171/Forsythe+-+A+model+comparison+for+daylength+as+a+function+of+latitude+and+day+of+year+-+1995.pdf
    # p=18 for astronomical twilight
    day = month/12*365.25+79
    theta = 0.2163108+2.*np.arctan(0.9671396*np.tan(0.00860*(day-186)))
    phi = np.arcsin(0.39795*np.cos(theta))
    arccosarg = (np.sin(p*np.pi/180.)+np.sin(latitude/180.*np.pi)*np.sin(phi))/(np.cos(latitude/180.*np.pi)*np.cos(phi))
    if abs(arccosarg)>=1.:
        return 0.0
    return 24./np.pi * np.arccos(arccosarg)

def get_stereographic_data(sims, latitude=0., month=0., hour=0., albedo=0.2, area=4., airmassCoeff=0.2, randomCoeff=0.5, elevation_cut = 0):
    # latitude in degrees
    # month in months from spring euquinox
    # hours in hours since midnight
    latitude = latitude/180.*np.pi 
    tilt = 23.4*np.sin(month/6.*np.pi)/180.*np.pi
    hour = hour/12.*np.pi
    xy, mag = [], []     
    for name in sims:
        sim = sims[name]
        sun = np.array([-1.4959787e+11,0,0]) # in m
        sun = rotY(sun, tilt)
        sun_n = sun/np.linalg.norm(sun)

        obs = np.array([REarth, 0, 0])
        obs = rotY(obs, -latitude)
        obs = rotZ(obs, hour)
        obs_n = obs/np.linalg.norm(obs)

        xyz = np.zeros((sim.N,3),dtype="float64")
        sim.serialize_particle_data(xyz=xyz)
        xyz = xyz[1:] # remove earth


        lit = np.linalg.norm(np.cross(xyz,sun_n),axis=1)>REarth

        xyz = xyz[lit]

        xyz_n = xyz/np.linalg.norm(xyz,axis=1)[:,np.newaxis]
        xyz_r = xyz - obs
        xyz_rd = np.linalg.norm(xyz_r,axis=1)
        xyz_rn = xyz_r/xyz_rd[:,np.newaxis]

        phase = np.arccos(np.clip(np.dot(xyz_rn, -sun_n), -1.0, 1.0)) # assume sun is in -x direction

        fac1 = 2/(3*np.pi**2)
        pfac = 3.1
        m_sun = -26.47 # g' band 
        magV = m_sun -2.5*np.log10(fac1 * area * albedo * ( (np.pi-phase)*np.cos(phase) + np.sin(phase) ) ) + 5 * np.log10(xyz_rd)
        #magV = m_sun -2.5*np.log10(2/(3*np.pi**(pfac+1)) * area * albedo * ( (np.pi-phase)*np.cos(phase) + np.sin(phase) )**pfac ) + 5 * np.log10(xyz_rd)


        elevation = (np.pi/2.-np.arccos(np.dot(xyz_rn,obs_n)))/np.pi*180.

        xyz = rotZ(xyz, -hour)
        xyz = rotY(xyz, latitude)
        xyz_r = xyz - np.array([REarth, 0, 0])
        xyz_rd = np.linalg.norm(xyz_r,axis=1)
        xyz_rn = xyz_r/xyz_rd[:,np.newaxis]

        #elevation_cut = 45
        xyz_rn = xyz_rn[elevation>elevation_cut]
        magV = magV[elevation>elevation_cut]

        airmass = getAirmass((90.-elevation[elevation>elevation_cut])*np.pi/180.)
        magV += airmassCoeff*airmass
        if randomCoeff>0.:
            magV += randomCoeff*np.random.normal(0.,1.,size=len(magV))

        xy.append(xyz_rn[:,1:3]/(1.+xyz_rn[:,0,np.newaxis]))
        mag.append(magV)
    if len(xy)>0:
        return np.concatenate(xy), np.concatenate(mag) 
    else:
        return None, None


def get_simulations(constellations=None, use_cache=True):
    if constellations is None:
        constellations = constellations_all
    sims = {}
    for c in constellations.keys():
        sim = None
        if use_cache:
            filename = "mega_"+"".join(x for x in c if x.isalnum())+".bin"
            try:
                sim = rebound.Simulation(filename)
            except:
                # need to create simulation
                pass
        if sim is None:
            sim = rebound.Simulation()
            sim.G = 6.67430e-11
            sim.add(m=MEarth)
            sim.N_active = 1
            add_to_simulation(sim, constellations[c])
            if use_cache:
                sim.save_to_file(filename)
        sims[c] = sim
    return sims
