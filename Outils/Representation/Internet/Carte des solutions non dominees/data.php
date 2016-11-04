<?php
// Variables globales
$request=array();
$stations=array();
$solutionNonDominees=array();






// Fonction pour generer les objets stations à partir du fichier de donnees .dat
function createXStations()
{
    global $stations;
    $myfile = fopen("grgs.pos+eop.990417.v60_MODIF_tran.dat", "r") or die("Unable to open file!");
    $numero=0;
    while (!feof($myfile))
    {
        if (fgets($myfile)==" DOMES_des_stations\n")
        {
            $ligne=fgets($myfile);
            while ($ligne!=" SOLN_des_stations\n")
            {
                $split=explode(" ", $ligne);
                $DOMES=trim($split[1]);
                $stations[$DOMES]=array("numero"=>$numero,"longitude"=>0,"latitude"=>0);
                $numero++;
                $ligne=fgets($myfile);
            }
        }
    }
    fclose($myfile);
}





// Fonction pour associer les coordonnées géographiques aux stations à partir du fichier .SNX
function loadLonLat()
{
    global $stations;
    $myfile = fopen("file_rfl0kr3sakkbj3k67lda3tok65_1452177263.SNX", "r") or die("Unable to open file!");
    while (!feof($myfile))
    {
        if (fgets($myfile)=="*CODE PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_\n")
        {
            $ligne=fgets($myfile);
            while ($ligne!="-SITE/ID                                                                        \n")
            {
                $split=preg_split('/[\s]+/', $ligne);
                if ($_POST["donneesStation"]==2)
                {
                    return $split ;
                }
                $DOMES=$split[3];
                if ($split[5]>=0)
                {
                    $longitude=$split[5]+$split[6]/60+$split[7]/3600;
                }
                else
                {
                    $longitude=$split[5]-$split[6]/60-$split[7]/3600;
                }
                
                if ($split[8]>=0)
                {
                    $latitude=$split[8]+$split[9]/60+$split[10]/3600;
                }
                else
                {
                    $latitude=$split[8]-$split[9]/60-$split[10]/3600;
                }
                $stations[$DOMES]["longitude"]=$longitude;
                $stations[$DOMES]["latitude"]=$latitude;
                $ligne=fgets($myfile);
            }
        }
    }
    fclose($myfile);
}





// Fonction qui récupere les solutions non dominees à partir du fichier .solutionNonDominees
function solutionNonDominees()
{
    global $solutionNonDominees;
    $myfile = fopen("grgs.pos+eop.990417.v60_MODIF.solutionNonDominees", "r") or die("Unable to open file!");
    while (!feof($myfile))
    {
        $ligne=fgets($myfile);
        if ($ligne!="")
        {
            $split=explode(" ", $ligne);
            $solutionNonDominees[$split[0]]=array("sigmax"=>$split[1],"sigmay"=>$split[2],"sigmaz"=>explode("\n", $split[3])[0]);
        }
    }
    fclose($myfile);
}






// Execution selon la requete ajax
if ($_POST["donneesStation"]==2)
{
    createXStations();
    echo json_encode(loadLonLat());
}
else if ($_POST["donneesStation"]==1)
{
    createXStations();
    loadLonLat();
    solutionNonDominees();
    $request["stations"]=$stations;
    $request["solutionsNonDominees"]=$solutionNonDominees;
    echo json_encode($request);
}
?>
