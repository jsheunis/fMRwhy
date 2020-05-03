
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
      "rest_run-1": "Rest - Run 1",
      "motor_run-1": "Motor - Run 1",
      "emotion_run-1": "Emotion - Run 1",
      "rest_run-2": "Rest - Run 2",
      "motor_run-2": "Motor - Run 2",
      "emotion_run-2": "Emotion - Run 2"
    };

    var str0 = "/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-qc/sub-001/func/sub-001_task-";
    var str1 = "/Users/jheunis/Desktop/sample-data/NEUFEPME_data_BIDS/derivatives/fmrwhy-qc/sub-001/func/PhysIO_task-";
    var str2 = "/sub-001_task-";
    var str3 = "_physioQC_03.jpg";
    var imgArray = {
      "mean_img": "_space-individual_mean.png",
      "std_img": "_space-individual_std.png",
      "tsnr_img": "_space-individual_tsnr.png",
      "grayplot_ro_img": "_echo-2_desc-RO_grayplot.png",
      "grayplot_gso_img": "_echo-2_desc-GSO_grayplot.png",
      "grayplot_lmotor_img": "_echo-2_desc-leftMotor_grayplot.png",
      "grayplot_biamy_img": "_echo-2_desc-bilateralAmygdala_grayplot.png"
    };
    // Reset sources for images
    for (var key in imgArray) {
      console.log("key " + key + " has value " + imgArray[key]);
      $("#" + key).attr("src", str0 + $(this).val() + imgArray[key]);
    }
    $("#physqcplots_img").attr("src", str1 + $(this).val() + str2 + $(this).val() + str3);

    // Reset heading names
    $("#spatialqcplots").html( "Spatial QC Plots: " + selectArray[$(this).val()] );
    $("#tempqcplots").html("Temporal QC plots: " + selectArray[$(this).val()] );
    $("#physqcplots").html("PhysIO QC plots: " + selectArray[$(this).val()] );

  });
});