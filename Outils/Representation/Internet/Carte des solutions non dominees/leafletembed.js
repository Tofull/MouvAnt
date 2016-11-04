// Creation de l'objet map
var map = L.map('map').setView([0,0], 1);





// Creation de variables generales pour l'affichage des marqueurs 
var markers = new L.FeatureGroup();
var listeMarker = [];





// Chargement des tuiles google
L.tileLayer('http://{s}.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',{
    maxZoom: 18,
    subdomains:['mt0','mt1','mt2','mt3']
}).addTo(map);
    
	
    
	

// Fonction qui permet d'ajouter des marqueurs sur la carte
function addMarkers (datum,data) {
    var markers = new L.FeatureGroup();
    var compteur=0;
    var listeStationsPrises=[];
    
    // Recuperation des stations sélectionnées par l'algorithme (valeur 1 dans le datum)
    for (var i = 0, len = datum.length; i < len; i++) {
        var character=datum[i];
        if (character=='1')
        {
            listeStationsPrises.push(compteur);
        }
        compteur++;
    }
    
    /* Transform the array in a dict */
    // Cela permet de tester si une station est selectionnée ou non en une seule ligne
    var listeStationsPrisesDico={}
    for(var i=0; i<listeStationsPrises.length;i++)
        {
            listeStationsPrisesDico[listeStationsPrises[i]]="";
        } 
        
    
    for (var station in data)
    {
        // Teste si une station est selectionnée ou non
        if ((data[station].numero in listeStationsPrisesDico))
        {
            var marker = L.marker([data[station].latitude, data[station].longitude]).bindPopup("Numero : "+station+"<br> latitude : "+data[station].latitude+"<br> longitude : "+data[station].longitude)
            listeMarker.push(marker);
            markers.addLayer(marker);
        }
    }
    
    // Ajout à la carte
    markers.addTo(map);
}

 



// Récuperation des donnees lorsque la page est chargée 
window.onload = function (e) {
    // Preparation de la requete ajax
	var ajax = new XMLHttpRequest(); 
	ajax.open('POST', 'data.php', true); 
	ajax.setRequestHeader('Content-type', 'application/x-www-form-urlencoded'); 
    // Declenchement de l'evenement lorsque le serveur renvoie des données
	ajax.addEventListener('readystatechange', 
		function(e) 
		{ 
			if(ajax.readyState == 4 && ajax.status == 200) 
			{ 
                // Recuperation des donnees
                var data = JSON.parse(ajax.responseText);
                // Affichage des donnees
                defile(data);
			} 
		}
	); 
    
    // Demande de donnees
	var data = "donneesStation=1"; // envoi de la requête
	ajax.send(data);
};           





// Fonction pour afficher les solutions non dominées
function defile(data)
{ 
    //suppression des anciens marqueurs
    markerDelAgain();
    
	//tirage au sort d'un datum 
    var nombre=0;
    var borneMax=Object.keys(data.solutionsNonDominees).length;
    nombre=Math.floor((Math.random() * borneMax));
    var compteur=0;
    for (var solution in data.solutionsNonDominees)
    {
        if (compteur==nombre)
        {
            var datum=solution;
            break;
        }
        compteur++;
    }    
    
    //affichage des stations du datum
    addMarkers(datum,data.stations);
    
    //animation
	window.setTimeout(function(){ defile(data); }, 500); 
}





/*Going through these marker-items again removing them*/
function markerDelAgain() {
    for(i=0;i<listeMarker.length;i++) 
    {
        map.removeLayer(listeMarker[i]);
    }  
}
