<script>
  let bg = -1;

  function getBg() {
    fetch("https://diabetessimapi.herokuapp.com/")
      .then(d => d.json())
    //   .then(d => (rand = d));
      .then(function (d) {
		  console.log(d)
		  return bg = d
	  });
  }
  
  function loadShow() {
	  getBg()
	  renderChart()
  }


// Copied from internet!
  	import { onMount } from "svelte";
	import Chart from "chart.js";
	onMount(async () => {});
	function renderChart() {
		var ctx = document.getElementById("myChart").getContext("2d");
		var chart = new Chart(ctx, {
		type: "line",
		data: {
			labels: [...Array(bg.length).keys()],
			datasets: [
			{
				label: "Blood glucose level",
				backgroundColor: "rgb(255, 99, 132)",
				borderColor: "rgb(255, 99, 132)",
				fill: 'False',
				data: bg
			}
			]
		},
		options: {}
		});

  }

// End copied from internet!
</script>

<!-- <h1>Your number is {bg}!</h1> -->
<!-- <button on:click={getBg}>Get a random number</button> -->

<h1>Press button to load standard BG values</h1>
<button on:click={loadShow}>Load</button> 
<canvas id="myChart"></canvas> 

<style>
:global(body){
	max-width: 1200px;
}
</style>