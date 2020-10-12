bids_dataset = {
  'sub': 'sub-001',
  'sessions': ['1'],
  'tasks': ['rest'],
  'runs': ['1'],
  'anat_template_session': '',
  'template_task': 'rest',
  'template_session': '',
  'template_run': '',
  'template_echo': '',
  'has_sessions': false,
  'has_runs': false,
  'is_multiecho': false,
  'map_rois': false,
};

function setupReportStructure() { 
  // alert(bids_dataset['sub']);
  // $('.anatroi').hide();

  // Get a reference to the table
  var tableRef = document.getElementById('kaaskoek');
  // Insert a row at the end of the table
  var newRow = tableRef.insertRow(-1);

  // Insert a cell in the row at index 0
  var newCell1 = newRow.insertCell(0);
  var newCell2 = newRow.insertCell(0);
  // Append a text node to the cell
  var newText1 = document.createTextNode('New category');
  var newText2 = document.createTextNode('New values');
  newCell1.appendChild(newText1);
  newCell2.appendChild(newText2);

  // If map_rois==false, remove nodes/elements related to ROIs
  if (!bids_dataset['map_rois']) {
    var element = document.getElementById('anatroitext');
    element.remove();
    var element2 = document.getElementById('anatroiimgs');
    // var node= document.getElementById("parent");
    element2.querySelectorAll('*').forEach(n => n.remove());
    element2.remove();
  }

  if (bids_dataset['has_sessions']) {

  } else {
    
  }
  
}


// $(document).ready(function () {
//   function createDOMelements(bids_dataset) {
//     alert("I am an alert box!");
//     if (!bids_dataset.sessions.isArray(array) || !bids_dataset.sessions.length) { // array does not exist, is not an array, or is empty
      
//     }

//     if (!bids_dataset.sessions.isArray(array) || !bids_dataset.sessions.length) { // array does not exist, is not an array, or is empty
      
//     }


//   }
// });



// document.addEventListener("DOMContentLoaded", function(event) { 
//   alert('Inside event listener in JS');
//   showHideROIs();
// });


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
// $(document).ready(function () {
//   $('#get-data').click(function () {
//     var showData = $('#show-data');

//     $.getJSON('assets/example.json', function (data) {
//       console.log(data);

//       var items = data.items.map(function (item) {
//         return item.key + ': ' + item.value;
//       });

//       showData.empty();

//       if (items.length) {
//         var content = '<li>' + items.join('</li><li>') + '</li>';
//         var list = $('<ul />').html(content);
//         showData.append(list);
//       }
//     });

//     showData.text('Loading the JSON file.');
//   });
// });


// Function to show/hide a set of divs based on dropdown selection
// $(document).ready(function() {
//   $('#colorselector').change(function(){
//     $('.colors').hide();
//     $('#' + $(this).val()).show();
//   });
// });


$(document).ready(function() {
  $('#runselector').change(function(){

    alert(bids_dataset['map_rois']);

    // Variables
    var selectArray = {
      "rest_run-1": "Rest - Run 1",
      "motor_run-1": "Motor - Run 1",
      "emotion_run-1": "Emotion - Run 1",
      "rest_run-2": "Rest - Run 2",
      "motor_run-2": "Motor - Run 2",
      "emotion_run-2": "Emotion - Run 2"
    };


    var str2 = "param_str2";
    var str3 = "param_str3";
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
      $("#" + key).attr("src", "img" + str2 + $(this).val() + imgArray[key]);
    }
    $("#physqcplots_img").attr("src", "img" + str2 + $(this).val() + str3);
    // Reset heading names
    $("#spatialqcplots").html( "Spatial QC Plots: " + selectArray[$(this).val()] );
    $("#tempqcplots").html("Temporal QC plots: " + selectArray[$(this).val()] );
    $("#physqcplots").html("PhysIO QC plots: " + selectArray[$(this).val()] );

  });
});