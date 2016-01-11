(function() {
  angular.module("Sunburst", [])
    .directive("sunburst", sunburst)
    .directive("onReadFile", onReadFile)
    .controller("MainCtrl", MainCtrl);

  // controller function MainCtrl
  function MainCtrl($http) {
    var ctrl = this;
    init();


    // function init
    function init() {
      // initialize controller variables
      ctrl.examples = [
        "out"
      ];
      ctrl.exampleSelected = ctrl.examples[0];
      ctrl.getData = getData;
      ctrl.selectExample = selectExample;
      
      // initialize controller functions
      ctrl.selectExample(ctrl.exampleSelected);
    }

    // function getData
    function getData($fileContent) {
      ctrl.data = $fileContent;
    }

    // function selectExample
    function selectExample(item) {
      var file = item + ".tab";
      $http.get(file).success(function(data) {
        ctrl.data = data;
      });
    }
  }


  // directive function sunburst
  function sunburst() {
    return {
      restrict: "E",
      scope: {
        data: "=",
      },
      link: sunburstDraw
    };
  }


  // directive function onReadFile
  function onReadFile($parse) {
    return {
      restrict: "A",
      scope: false,
      link: function(scope, element, attrs) {
        var fn = $parse(attrs.onReadFile);
        element.on("change", function(onChangeEvent) {
          var reader = new FileReader();
          reader.onload = function(onLoadEvent) {
            scope.$apply(function() {
              fn(scope, {
                $fileContent: onLoadEvent.target.result
              });
            });
          };
          reader.readAsText((onChangeEvent.srcElement || onChangeEvent.target).files[0]);
        });
      }
    };
  }
})();
