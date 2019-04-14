function Data=datainput_Schoen_17_2_5
%17-dimensional Schoen function
%--------------------------------------------------------------------------
%Copyright (c) 2012 by Juliane Mueller
%
% This file is part of the surrogate model module toolbox.
%
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%Tampere University of Technology, Finland
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%

Data.xlow=zeros(1,17); %lower variable bounds
Data.xup=ones(1,17); %upper variable bounds
Data.objfunction=@(x)myfun(x); %handle to objective function
Data.dim=17;
Data.integer=[];
Data.continuous=(1:17);

end %function

function y = myfun(x)
x=x(:)';

StatPoint = [...
0.86386309049428855000  0.50123228440306045000  0.71169251617075890000  0.86095369523449139000  0.38723968193363201000  0.77473309497706078000  0.36465071812102612000  0.19886784709287861000  0.74419563357013241000  0.21149041120727335000  0.82486249284327251000  0.81335213025173925000  0.12340990026968450000  0.40684364351772451000  0.38497729147412107000  0.83575312899054621000  0.63423568689491450000  
0.90618764229090443000  0.50469551299110338000  0.78368123702386527000  0.91689722280854602000  0.77521797931656133000  0.74010261123574461000  0.64745801459331265000  0.80310383821559905000  0.32510055638289648000  0.33514742904836359000  0.27314548957826384000  0.97827786159267804000  0.77250702096139123000  0.90168923040767934000  0.92060008725583031000  0.29921109887401209000  0.71527665222725223000  
0.20014825763326596000  0.98376455787427119000  0.38713417299738945000  0.66584966126015810000  0.32382761435164076000  0.41654891645667930000  0.72693100731370952000  0.09832448442242547300  0.16743688173199231000  0.75225601758075367000  0.64123586202875082000  0.28658334339100960000  0.95037304635058339000  0.96303205251622381000  0.17684661995686121000  0.04128515618399991600  0.46598936434952903000  
0.85939524462774308000  0.84229176307863374000  0.46375642685563229000  0.39694172223292118000  0.04354313025720402100  0.30178370216390515000  0.07936138855480748300  0.47662571397337572000  0.38311685884624780000  0.24500089576542450000  0.97556765454907568000  0.94669318633800070000  0.64864845648517899000  0.50591317999422125000  0.13027972153564008000  0.39197572959923971000  0.64840008929870441000  
0.79595557163994901000  0.45119253993011704000  0.92187252769159400000  0.32962646281355040000  0.90835668382678647000  0.43536612036944322000  0.35556992389456910000  0.32274706515807355000  0.87073185730277847000  0.15996848278403078000  0.78261725968681251000  0.45656350426587439000  0.10139689451975499000  0.62819917542031389000  0.13802729981505124000  0.62295196969578115000  0.31902653299303735000  
0.33819772500615269000  0.31997624506494488000  0.50442374129719336000  0.55114667448764665000  0.87508521922178983000  0.08339342789688375400  0.47542585196294002000  0.49580080906358187000  0.43940398665982938000  0.94999329399623555000  0.08786461212417914800  0.42776089674367868000  0.42975917946859166000  0.51453541584634754000  0.44290089148830653000  0.78099413937985540000  0.20940317791395724000  
0.72932294081657323000  0.23444362544409550000  0.05305357422825739100  0.60039602077558929000  0.70388541560791307000  0.56491842258351188000  0.08150361692147963000  0.38688149538226685000  0.73870639879930533000  0.00525882123700986690  0.70361081388712410000  0.56625381924250129000  0.12382495471624808000  0.88399584918297702000  0.58254706383108457000  0.21606016592136720000  0.43202950531189044000  
0.54275043625016162000  0.76816047543892663000  0.77176863559443576000  0.06615864364561166400  0.97525644771666342000  0.30980059483802741000  0.38032997752226949000  0.78055555307832836000  0.58297418638498544000  0.05290397364662984700  0.04091483359181833500  0.02992391890265526700  0.98482746894903594000  0.32655129752041423000  0.97773421215741574000  0.44007567902987599000  0.24725822692525048000  
0.43020085578946105000  0.32136705846610392000  0.84376409798316310000  0.10664142386500625000  0.26690434141061831000  0.00373731958623049410  0.37261211808392114000  0.78976570211198527000  0.77563866383438718000  0.92759773342123142000  0.66486667518948928000  0.49456241532662115000  0.75837873493405450000  0.48905274507098628000  0.99637655863501773000  0.00135469946095211850  0.49071216950696694000  
0.39633835868284079000  0.55326268794466749000  0.70546632825036959000  0.21181755084433523000  0.67633173927893298000  0.08357596126434203900  0.25780431570301382000  0.46236056214213028000  0.24595152961973143000  0.67863624593924987000  0.79979532680922039000  0.24158542287849874000  0.27002186516398952000  0.41847502539423259000  0.21545465131263478000  0.51566730171936748000  0.15142580968069383000  
0.21374603030925002000  0.91682553892766505000  0.20784672765714854000  0.27212960115131579000  0.21354674179320721000  0.70234082315803714000  0.16656009353657653000  0.40320670704681522000  0.93849130809129233000  0.38990914470537524000  0.86065969785971996000  0.83315793250277537000  0.92973899307830632000  0.64046602638126693000  0.40744163813536982000  0.30737442082209659000  0.56935026611133532000  
0.13375768119377171000  0.26490840823402084000  0.26833997232956502000  0.79261866520294821000  0.15777721897542724000  0.86610732183797690000  0.24991397714989355000  0.69948564275204361000  0.06853889718898136000  0.27996407205089602000  0.76287160435426227000  0.78538454464600338000  0.43398347734793963000  0.12268995194395471000  0.60573620849240084000  0.79415491549090611000  0.96292900639541135000  
0.88272295072251228000  0.47875877103877729000  0.72148595890392475000  0.26431903821954456000  0.94544859093191946000  0.86270186048287789000  0.15355821905176689000  0.55107091124254548000  0.29680431746130981000  0.64499206073820470000  0.20467983150731264000  0.69201176837948075000  0.67174389951792934000  0.06157637399918425300  0.55913176589784819000  0.45462800983644186000  0.16731745478491183000  
0.81474822695928761000  0.51165546021301656000  0.74204872109574194000  0.46511536679706411000  0.82534248947782074000  0.04434047645023940500  0.61554192117382833000  0.59003073463187494000  0.70493881093070254000  0.51764943229821880000  0.34963112694264969000  0.02788210467971363400  0.89989741099533116000  0.75175840970040719000  0.10505821640512374000  0.89194336637580263000  0.69427801569948455000  
0.74534506606984929000  0.94917934428612794000  0.85724220928972983000  0.01053834261565120600  0.36025060033102879000  0.71220037653852819000  0.04407436732216109300  0.77695010345233217000  0.27850136262939534000  0.02780962060327596800  0.50271757456099875000  0.25833714023180931000  0.78027050519085328000  0.52374295229735335000  0.51629612405490766000  0.14246454175825998000  0.03963951395951104800  
0.97604521850329529000  0.32638953865785070000  0.53632471054654129000  0.61933187724165950000  0.51931972622173528000  0.11722378548849502000  0.06837125900887014400  0.04088064622858421400  0.66053466005181649000  0.16968716264559450000  0.02767362452183503400  0.90422577687001038000  0.34659339762135355000  0.25392697329773500000  0.89542615097552314000  0.97796390387625332000  0.96036534597799572000  
0.12427353951456570000  0.05931050421674607900  0.78542758872859708000  0.54167943277697439000  0.90316688396574696000  0.72290410659221316000  0.30422134027171649000  0.61382527439436640000  0.36134029877827345000  0.01101057080740435600  0.68471921612462572000  0.26598153235230820000  0.74489686309862457000  0.40479851999258487000  0.68042687659755474000  0.90424002936228087000  0.96776849967824241000  
0.63125497340963310000  0.13085266294875428000  0.67947932547490342000  0.45388815983244124000  0.99385816493409773000  0.68621249884014934000  0.56451424314784671000  0.12590427979459692000  0.12343773865517797000  0.98304964527425154000  0.99014041752220971000  0.99902874912758466000  0.84721915307435280000  0.63365012888838212000  0.75547725066057425000  0.22862007253030589000  0.71263413035071543000  
0.93099425239086153000  0.59373646237389721000  0.81719432037718054000  0.71177932772674524000  0.78148677687833235000  0.94252378507445489000  0.32424004225313308000  0.24254899995519347000  0.82659030916040532000  0.98365076240340288000  0.83998633572162307000  0.16557218399277171000  0.48854507589378621000  0.39016652804165736000  0.75250414039966229000  0.16780644910814238000  0.71490277568358029000  
0.86551811359721120000  0.88007918651998662000  0.39119249492616592000  0.59731865624807812000  0.28018444089027866000  0.63946690404041528000  0.48417942349920529000  0.36396979717030353000  0.06496936304152509100  0.59473085930017988000  0.36596117176645165000  0.82730998368842712000  0.18429748660621825000  0.72269712002233422000  0.61796036803021703000  0.51115686412611705000  0.73712356338505902000];

StatValue = [...
100.00000000000000000000
42.87585077488715500000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
23.01398420133861200000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000
100.00000000000000000000];

k = length(StatValue);
Xtemp = ones(k,1)*x;
TempNorm = sum((Xtemp-StatPoint).^2,2);
TempProd = ones(1,k);
for i = 1:k
    TempProd(i) = prod(TempNorm(1:(i-1)))*prod(TempNorm((i+1):k));
end
y = (TempProd*StatValue)/sum(TempProd);

end %myfun