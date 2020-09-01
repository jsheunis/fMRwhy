
function openTab1(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent1");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks1");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}


function openTab2(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent2");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks2");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}

function openTab3(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent3");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks3");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}

function openTab4(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent4");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks4");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}


function openTab5(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent5");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks5");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}


function openTab6(evt, tabName) {
  // Declare all variables
  var i, tabcontent, tablinks;

  // Get all elements with class="tabcontent" and hide them
  tabcontent = document.getElementsByClassName("tabcontent6");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }

  // Get all elements with class="tablinks" and remove the class "active"
  tablinks = document.getElementsByClassName("tablinks6");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }

  // Show the current tab, and add an "active" class to the button that opened the tab
  document.getElementById(tabName).style.display = "block";
  evt.currentTarget.className += " active";
}






// function to load data from external json file (does not work in chrome since Chrome does not have acces to local filesystem)
$(document).ready(function () {
  $('#get-data').click(function () {
    var showData = $('#show-data');

    $.getJSON('assets/example.json', function (data) {
      console.log(data);

      var items = data.items.map(function (item) {
        return item.key + ': ' + item.value;
      });

      showData.empty();

      if (items.length) {
        var content = '<li>' + items.join('</li><li>') + '</li>';
        var list = $('<ul />').html(content);
        showData.append(list);
      }
    });

    showData.text('Loading the JSON file.');
  });
});


// Function to show/hide a set of divs based on dropdown selection
$(document).ready(function() {
  $('#colorselector').change(function(){
    $('.colors').hide();
    $('#' + $(this).val()).show();
  });
});


$(document).ready(function() {
  $('#runselector').change(function(){

    // Variables
    var selectArray = {
      "motor_run-1": "Motor - Run 1",
      "emotion_run-1": "Emotion - Run 1",
      "rest_run-2": "Rest - Run 2",
      "motor_run-2": "Motor - Run 2",
      "emotion_run-2": "Emotion - Run 2"
    };

    var statsArray = {
      "motor_run-1": "FingerTapping",
      "emotion_run-1": "Faces>Shapes",
      "motor_run-2": "MentalFingerTapping",
      "emotion_run-2": "MentalEmotion"
    };


    var str2 = "/sub-001_task-";
    var str3 = "_physioQC_03.jpg";
    var imgArray = {
      "bold1_img": "_echo-1_desc-rapreproc_bold.png",
      "bold2_img": "_echo-2_desc-rapreproc_bold.png",
      "bold3_img": "_echo-3_desc-rapreproc_bold.png",
      "bold4_img": "_desc-combinedMEtsnr_bold.png",
      "bold5_img": "_desc-combinedMEt2star_bold.png",
      "bold6_img": "_desc-combinedMEte_bold.png",
      "bold7_img": "_desc-combinedMEt2starFIT_bold.png",
      "bold8_img": "_desc-t2starFIT_bold.png",
      "tsnr1_img": "_echo-2_desc-rapreproc_tsnr.png",
      "tsnr2_img": "_desc-combinedMEtsnr_tsnr.png",
      "tsnr3_img": "_desc-combinedMEt2star_tsnr.png",
      "tsnr4_img": "_desc-combinedMEte_tsnr.png",
      "tsnr6_img": "_desc-combinedMEt2starFIT_tsnr.png",
      "tsnr7_img": "_desc-t2starFIT_tsnr.png",
      "tsnr5_img": "_desc-tsnrPercdiffRainclouds_brain.png",
      "percdiff1_img": "_desc-combinedMEtsnr_percdiff.png",
      "percdiff2_img": "_desc-combinedMEt2star_percdiff.png",
      "percdiff3_img": "_desc-combinedMEte_percdiff.png",
      "percdiff4_img": "_desc-combinedMEt2starFIT_percdiff.png",
      "percdiff5_img": "_desc-t2starFIT_percdiff.png",
      "roi1_img": "_echo-2_desc-srapreproc_tsplot.png",
      "roi2_img": "_desc-scombinedMEtsnr_tsplot.png",
      "roi3_img": "_desc-scombinedMEt2star_tsplot.png",
      "roi4_img": "_desc-scombinedMEte_tsplot.png",
      "roi6_img": "_desc-scombinedMEt2starFIT_tsplot.png",
      "roi7_img": "_desc-st2starFIT_tsplot.png",
      "roi5_img": "_desc-tsnrPercdiffRainclouds_roi.png",
      "stats1_img": "_echo-2_desc-" + statsArray[$(this).val()] + "_threshtmap.png",
      "stats2_img": "_echo-combinedMEtsnr_desc-" + statsArray[$(this).val()] + "_threshtmap.png",
      "stats3_img": "_echo-combinedMEt2star_desc-" + statsArray[$(this).val()] + "_threshtmap.png",
      "stats4_img": "_echo-combinedMEte_desc-" + statsArray[$(this).val()] + "_threshtmap.png",
      "stats9_img": "_echo-combinedMEt2starFIT_desc-" + statsArray[$(this).val()] + "_threshtmap.png",
      "stats10_img": "_echo-t2starFIT_desc-" + statsArray[$(this).val()] + "_threshtmap.png",
      "stats5_img": "_echo-2_desc-" + statsArray[$(this).val()] + "_noFWEthreshtmap.png",
      "stats6_img": "_echo-combinedMEtsnr_desc-" + statsArray[$(this).val()] + "_noFWEthreshtmap.png",
      "stats7_img": "_echo-combinedMEt2star_desc-" + statsArray[$(this).val()] + "_noFWEthreshtmap.png",
      "stats8_img": "_echo-combinedMEte_desc-" + statsArray[$(this).val()] + "_noFWEthreshtmap.png",
      "stats11_img": "_echo-combinedMEt2starFIT_desc-" + statsArray[$(this).val()] + "_noFWEthreshtmap.png",
      "stats12_img": "_echo-t2starFIT_desc-" + statsArray[$(this).val()] + "_noFWEthreshtmap.png",
    };

    // Reset sources for images
    for (var key in imgArray) {
      console.log("key " + key + " has value " + imgArray[key]);
      $("#" + key).attr("src", "img" + str2 + $(this).val() + imgArray[key]);
    }
    $("#physqcplots_img").attr("src", "img" + str2 + $(this).val() + str3);
    // Reset heading names
    $("#boldtabsheading").html( "A. BOLD: " + selectArray[$(this).val()] );
    $("#tsnrtabsheading").html( "B. tSNR: " + selectArray[$(this).val()] );
    $("#percdifftabsheading").html( "C. Percentage difference in tSNR: " + selectArray[$(this).val()] );
    $("#roitabsheading").html( "D. Region of interest: " + selectArray[$(this).val()] );
    $("#statstabsheading").html( "E. Statistical maps: " + selectArray[$(this).val()] );



  });
});
