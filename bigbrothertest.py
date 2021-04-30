import astroplan
import astropy
import numpy
import os
import requests
from lxml import html
import re
import datetime

import datetime
import astropy.coordinates
import astropy.time
import astropy.units as u
from astropy.time import TimezoneInfo
import ephem
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter, attrgetter
from random import randint
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from rts2solib import rts2comm, queue
from rts2solib import scriptcomm

#what will happen:
#	a trigger will come in
#	this trigger will have some sort of
#		localization: set of RA and DEC
#					  priorities
global colors
colors = [
	"#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
	"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
	"#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
	"#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
	"#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
	"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
	"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
	"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
	"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
	"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
	"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
	"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
	"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
	"#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
	"#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
	"#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"
	]

class trigger:
	def __init__(self):
		self.trigger_requests = []
		self.date = None
		self.source = None


class trigger_request:
	def __init__(self,name, ra, dec, priority):
		self.ra = ra
		self.dec = dec
		self.priority = priority
		pass


class observatory:
	def __init__(self):
		self.name = None
		self.latitude = None
		self.longitude = None
		self.status = None
		self.telescopes = []
		pass


class telescope:
	def __init__(self):
		self.name = None
		self.elevation = None
		self.instrument = None
		self.status = None
		self.queue = []
		pass


#class queue:
#	def __init__(self):
#		self.date = None
#		self.queue_objs = []
#		pass


class queue_obj:
	def __init__(self, name, ra, dec, priority, exp_objs, constraints=None, rts2ids=[]):
		self.name = name
		self.ra = ra
		self.dec = dec
		self.priority = priority
		self.exp_objs = exp_objs
		self.rts2ids = rts2ids
		#self.skycoord = astropy.coordinates.SkyCoord(
		#	ra=astropy.coordinates.Angle(self.ra),
		#	dec=astropy.coordinates.Angle(self.dec)
		#)
		self.skycoord = astropy.coordinates.SkyCoord(
			'{} {}'.format(self.ra, self.dec),
			unit=(u.hourangle, u.deg)
		)
		self.overhead = 0
		self.estimateoverhead()

		self.constraints = constraints
		self.constrained_visible = False

		self.altaz = None
		self.frame_airmass = None
		self.peak_airmass = None

		self.start_observation = None
		self.end_observation = None
		self.scheduled = False

	def set_altaz(self, frame, setpeak=False):
		self.altaz = self.skycoord.transform_to(frame)
		airmasses = []
		for aa in self.altaz.secz:
			if aa > 1:
				airmasses.append(aa)
		print(self.name, len(airmasses))
		if len(airmasses):
			if setpeak:
				self.peak_airmass = min(airmasses)

			self.frame_airmass = min(airmasses)
			self.airmass_potential = self.peak_airmass-self.frame_airmass
		else:
			self.peak_airmass = 99
			self.frame_airmass = 99
			self.airmass_potential = 99


	def outputobjectinfo(self):
		print("Queue Object: {}".format(self.name))
		print("RA: {}".format(self.ra))
		print("DEC: {}".format(self.dec))
		print("Est. Overhead: {}".format(self.overhead))
		print("Observation Infos")
		for obsinfo in self.exp_objs:
			print("     Filter: {}, Exposure Time: {}, Amount {}".format(obsinfo.filter, obsinfo.exptime, obsinfo.amount))

	def estimateoverhead(self, readout=10, slew_rate=60):
		#values are in seconds
		#initial overhead is the slewrate: 60 seconds
		t = slew_rate

		#add up all of the exposure times: sum(exposuretime*amount)
		exp_time = sum(float(x.exptime)*float(x.amount) for x in self.exp_objs)

		#apply readout times for every single exposure: readouttime*N_exposure
		read_outs = readout*sum(int(x.amount) for x in self.exp_objs)

		#add it all together
		t = t+exp_time+read_outs

		self.overhead = t
		return t

	def evaluate_constraints(self, time, location):
		if self.constraints is not None:
			if self.constraints.moon_distance_threshold is not None:
				moon = astropy.coordinates.get_moon(time=time, location=location)
				sep = self.skycoord.separation(moon)/u.deg
				#print('Evaluating constraints moon distance: {}, {}>{}'.format(
				#	self.name,
				#	sep,
				#	self.constraints.moon_distance_threshold
				#))
				if sep > self.constraints.moon_distance_threshold:
					self.constrained_visible = True
				else:
					print('{} moon distance constrained! {}'.format(self.name, sep))
					self.constrained_visible = False
					return
			if self.constraints.airmass_threshold is not None:
				#print('Evaluating constraints airmass: {}, {}<{}'.format(
				#	self.name,
				#	self.frame_airmass,
				#	self.constraints.airmass_threshold
				#))
				if self.frame_airmass < self.constraints.airmass_threshold:
					self.constrained_visible = True
				else:
					print('{} airmass constrained! {}'.format(self.name, self.frame_airmass))
					self.constrained_visible = False
					return
		else:
			#no constraints -> it is visible
			self.constrained_visible = True


class exp_obj:
	def __init__(self, amount, filter, exptime):
		self.amount = amount
		self.filter = filter
		self.exptime = exptime


class queue_obj_constraints:
	def __init__(self, moon_distance_threshold=10, airmass_threshold=2.0):
		self.moon_distance_threshold = moon_distance_threshold #arcseconds?
		self.airmass_threshold = airmass_threshold

def findoffset(sr):
	type_dict = {"UVOT":0, "AzTEC":1, "SPOL":1, "STAND":3}
	type_i = 14
	for key in type_dict.keys():
		if key in sr:
			typeindex = sr.index(key)
			if key == "STAND":
				return -99
			return type_i - typeindex
	return 0


def formatcoord(coord):
	sine=""
	if "-" in coord:
		coord = coord.split("-")[1]
		return "-{}d{}m{}s".format(coord[0:2],coord[2:4],coord[4:])
	return  "{}d{}m{}s".format(coord[0:2],coord[2:4],coord[4:])



def orp_session():
    orpdbpath = 'postgresql+psycopg2://artn:ArTn_520@scopenet.as.arizona.edu'
    engine = create_engine(orpdbpath)
    meta = MetaData()
    meta.reflect(bind=engine)
    session = sessionmaker(bind=engine)()

    return meta, session


def orp_targets(queued_iso=datetime.datetime(2021, 3, 22)):
    
    queued_iso_end = queued_iso + datetime.timedelta(days=6)
    meta, session = orp_session()
    obsreqs = meta.tables['obsreqs']
    obs_col = meta.tables['obsreqs'].columns

    db_resp = session.query(obsreqs).filter(
            obs_col['queued']==True,
            #obs_col['completed']==False,
            obs_col['queued_iso'] > queued_iso,
            obs_col['queued_iso'] < queued_iso_end
        ).all()

    '''Dirty way to group obsreqs'''
    target_names = []
    for resp in db_resp:
        objname = resp.object_name
        if '.us.' in objname:
            objname = objname.split('.us.')[0]
        target_names.append(objname)
    target_names = list(set(target_names))

    return_data = []
    for t in target_names:
        observations = [x for x in db_resp if t in x.object_name]
        ob1 = observations[0]
        observation_infos = []
        rts2ids = []
        for ob in observations:
            rts2ids.append(int(ob.rts2_id))
            observation_infos.append(
                    exp_obj(
                        ob.num_exp,
                        ob.filter_name,
                        ob.exp_time
                    )
                )
        q = queue_obj(
                t,
                ob1.ra_hms,
                ob1.dec_dms,
                5,
                observation_infos,
                queue_obj_constraints(
                    moon_distance_threshold=10, 
                    airmass_threshold=float(ob1.airmass)    
                ),
                rts2ids = rts2ids
            )
        return_data.append(q)

    return return_data


def readfromweb():
	page = requests.get("http://slotis.kpno.noao.edu/LOTIS/skypatrol.html")
	tree = html.fromstring(page.content)
	targetsinfo = tree.xpath('//h2/text()')[3:]
	targetsinfo = [re.sub('\n', '', x) for x in targetsinfo if '%' in x]
	lotisdata = []
	for line in targetsinfo:
		splitline = line.split()
		linelen = len(splitline)
		os = findoffset(splitline)
		if os != -99:
			name = splitline[9]
			ra = formatcoord(splitline[2])
			dec = formatcoord(splitline[3])
			exptime = splitline[5]
			filters = splitline[7]
			amounts = splitline[6]
			objtype = ""
			priority = randint(1,5)
			observationinfos = []
			for f in filters:
				observationinfos.append(exp_obj(amounts, f, exptime))
			lotisdata.append(queue_obj(name,ra,dec,priority,observationinfos,queue_obj_constraints()))
	return lotisdata

def aztecTargs03222021():
	targs = [
		queue_obj(
			'SN2021epp',
			'08:10:55.29', '-06:02:49.4',
			5,
			[
				exp_obj(4, 'V', 200),
				exp_obj(4, 'U', 300),
				exp_obj(4, 'B', 300),
				exp_obj(4, 'R', 200),
				exp_obj(4, 'I', 200)
			],
			queue_obj_constraints()
		),
		queue_obj(
			'SN2021dov',
			'08:56:18.71','-00:26:31.7',
			5,
			[
				exp_obj(4, 'U', 90),
				exp_obj(4, 'B', 90),
				exp_obj(4, 'V', 60),
				exp_obj(4, 'R', 45),
				exp_obj(4, 'I', 60)
			],
			queue_obj_constraints()
		),
		queue_obj(
			'SN2020adow',
			'08:33:42.18','+27:42:43.3',
			5,
			[
				exp_obj(4, 'U', 200),
				exp_obj(4, 'B', 200),
				exp_obj(4, 'V', 180),
				exp_obj(4, 'R', 120),
				exp_obj(4, 'I', 120)
			],
			queue_obj_constraints()
		),
		queue_obj(
			'SN2021dbg',
			'09:24:26.80','-06:34:53.1',
			5,
			[
				exp_obj(4, 'U', 120),
				exp_obj(4, 'B', 120),
				exp_obj(4, 'V', 90),
				exp_obj(4, 'R', 90),
				exp_obj(4, 'I', 90)
			],
			queue_obj_constraints()
		),
		queue_obj(
			'SN2021bnw',
			'10:53:52.18','+12:33:29.1',
			5,
			[
				exp_obj(4, 'U', 200),
				exp_obj(4, 'B', 200),
				exp_obj(4, 'V', 180),
				exp_obj(4, 'R', 120),
				exp_obj(4, 'I', 120)
			],
			queue_obj_constraints()
		),
		queue_obj(
			'SN2021gmj',
			'10:38:47.17','+53:30:31.0',
			5,
			[
				exp_obj(5, 'U', 240),
				exp_obj(5, 'B', 240),
				exp_obj(5, 'V', 180),
				exp_obj(5, 'R', 150),
				exp_obj(5, 'I', 180)
			],
			queue_obj_constraints()
		)
	]
	return targs

class FFTarget():

    def __init__(self, targ, frame):
        self.id = targ[0][0]
        self.name = targ[0][1]
        self.ra = targ[0][2]
        self.dec = targ[0][3]

        self._setSkyCoord()
        self._setAltAz(frame)

    def _setSkyCoord(self):
        self.skycoord = astropy.coordinates.SkyCoord(ra = self.ra*u.degree, dec = self.dec*u.degree)

    def _setAltAz(self, frame):
        self.altaz = self.skycoord.transform_to(frame)

    def dump(self):
        for key, value in self.__dict__.items():
            print('{}: {}, {}'.format(key, value, type(value)))

def FocusRunDecider(frame, frame_night):
    FOCUS_FIELD_IDS = [3611, 3612, 3613, 3614, 3615, 3616, 3617, 3618, 3619, 3620, 3621, 3622, 3623, 3624]
    com = rts2comm()
    focus_field_targets = {}

    for _id in FOCUS_FIELD_IDS:
        t = com.get_target(_id)
        fft = FFTarget(t, frame)
        secz = fft.altaz.secz

        if secz > 0.0 and secz < 3:
            print(secz, fft.name)
            focus_field_targets[str(secz)] = fft

    lowest_secz = min(list(focus_field_targets.keys()))
    focus_field = focus_field_targets[lowest_secz]
    
    return_field = queue_obj(
        focus_field.name,
        focus_field.ra, focus_field.dec,
        0,
        [
                exp_obj(2, 'V', 30),
        ],
        queue_obj_constraints(
            moon_distance_threshold=1,
            airmass_threshold=float(3.0)
        ),
        rts2ids = [focus_field.id]
    )
    return_field.skycoord = focus_field.skycoord
    return_field.set_altaz(frame_night, setpeak=True)
    
    return return_field


def observing_schedule(targets, nightrange, midnight, location, frame_night):

	#First method:
	#	split the observability night range into a number of hours
	#	for each range:
	#		find all objects that are observable
	#			eliminate those that are constrained
	#			sort them by priority
	#			then by highest airmass
	#		Once sorted. Start with the first
	#			Set start_observation = midnight - nightrange[0]*u.hr
	#			Set end_observation = self.start_observation + self.overhead
	#			Set queued = True

	print("Starting Scheduler!")

	#The range of hours of observability
	hour_intervals = int(nightrange[len(nightrange)-1]-nightrange[0])
	num_slices = int(len(nightrange)/hour_intervals)
	frames_range = nightrange*u.hour
	print(hour_intervals, num_slices)

	#set the first observation time to be the first timeslot 
	#	of astronomical night
	obs_time = midnight + frames_range[0]
	print(obs_time)
	
	unscheduled_targets = targets
	scheduled_targets = []

	for ni in range(hour_intervals):
		#get all unscheduled targets

		unscheduled_targets = [t for t in unscheduled_targets if not t.scheduled]
		hour = ni+1           

		#find the hour intervals
		#	if it is the last interval
		#		get the whats left after the last slice
		#	else: get the appropriate slice
		if hour == hour_intervals:
			frame_hours = midnight + frames_range[ni*num_slices:]
		else:
			frame_hours = midnight + frames_range[ni*num_slices:hour*num_slices]

		#create a new frame based on the frame_hour slice
		#	and evaluate the altaz for each target in this frame
		frame = astropy.coordinates.AltAz(
			obstime=frame_hours,
			location=location
		);
		focus_frame = astropy.coordinates.AltAz(
			obstime=frame_hours[0],
			location=location
		);

		if ni == 0 or ni == int(hour_intervals/2.0):
			focus_field = FocusRunDecider(focus_frame, frame_night)
			#focus_field.set_altaz(frame_night, setpeak=True)
			unscheduled_targets.append(focus_field)

		[t.set_altaz(frame) for t in unscheduled_targets]

		#do something with constraints: boolean observable = T/F?
		[t.evaluate_constraints(obs_time, location) for t in unscheduled_targets]
		constrained_targets = [t for t in unscheduled_targets if t.constrained_visible]

		#Scheduling logic:
		#	loop over priorities
		#	group "constrained_visible" targets in that priority
		#	sort those grouped targets by priority then by airmass
		#	set thos at the top to be scheduled
		#	-> do it again through all priorities

		priority_range = [0,1,2,3,4,5]

		for p in priority_range:

			if obs_time > frame_hours[-1]:
					break

			print('Priority Group: {}'.format(p))
			#This will group the targets based on priority
			#	allows for decimal priority as well
			priority_targs = [t for t in constrained_targets if t.priority > p-1 and t.priority <= p]

			#Sort the targets by their airmasss potential.
			sort_p_targs = sorted(priority_targs, key=attrgetter('airmass_potential'), reverse=True)
			for s in sort_p_targs:
				print("   {}: {}, {}, {}".format(s.name, s.frame_airmass, s.priority, s.constrained_visible))

			#loop over each of these sorted targets
			#	assigning the top target an observation slot
			for s in sort_p_targs:
				#set the start_time to be the current obs_time
				#	end_time to be the obstime+overhead
				#	set to be scheduled
				s.start_observation = obs_time
				s.end_observation = obs_time + s.overhead*u.second
				s.scheduled = True

				print('{} scheduled: {} with overhead: {} at Priority {}'.format(
					s.name,
					s.frame_airmass,
					#s.start_observation,
					#s.end_observation,
					s.overhead,
					s.priority
				))
				scheduled_targets.append(s)

				#reset the current obs_time + smol overhead
				obs_time = s.end_observation+(3*u.minute)

				#if the obs_time is greater than the last time in the frame:
				#	break out of the priority loop and create a new frame block
				if obs_time > frame_hours[-1]:
					print("Iterating new frame")
					break

	#give back the scheduled targets
	return(scheduled_targets)

def main():

	fig = plt.figure(figsize=(11,7))

	big2 = astropy.coordinates.EarthLocation.of_site('mtbigelow')

	targets = orp_targets()#queued_iso=datetime.datetime.now())
	#targets = aztecTargs03222021()targets = orp_targets()


	utcoffset = 7*astropy.units.hour

	d = datetime.datetime.now()
	#d = datetime.datetime(2021, 3, 22)

	midnight = astropy.time.Time('{}-{}-{} 00:00:00'.format(d.year, d.month, d.day)) + utcoffset
	
	delta_midnight = np.linspace(-12, 12, 1000)
	day_times = midnight + delta_midnight*u.hour

	frame_night = astropy.coordinates.AltAz(
		obstime=day_times,
		location=big2
	)

	#set the targets peak airmass:
	[t.set_altaz(frame_night, True) for t in targets]

	[print(t.peak_airmass) for t in targets]

	sun_altaz = astropy.coordinates.get_sun(day_times).transform_to(frame_night)
	moon_altaz = astropy.coordinates.get_moon(day_times).transform_to(frame_night)

	sun_alt_beforemid = []
	nightrange = []
	for i,x in enumerate(sun_altaz.alt < 5*u.deg):
		if x:
			nightrange.append(delta_midnight[i])
			sun_alt_beforemid.append(x)

	astronomical_night = []
	for i,x in enumerate(sun_altaz.alt < -18*u.deg):
		if x:
			astronomical_night.append(delta_midnight[i])

	scheduled_targets = observing_schedule(targets=targets, nightrange=astronomical_night, midnight=midnight, location=big2, frame_night=frame_night)

	moon_delta_midnight = []
	moon_airmass = []
	for i,x in enumerate(moon_altaz.secz):
		if x <= 5 and x > 1.0:
			moon_delta_midnight.append(delta_midnight[i])
			moon_airmass.append(x)

	#plotting the beginning and ending of the astronomical night
	plt.axvline(x=astronomical_night[0], color='k', linestyle='-', linewidth=2)
	plt.text(astronomical_night[0]-0.7, 0.93, 'Night begin', fontsize=10)
	plt.axvline(x=astronomical_night[-1], color='k', linestyle='-', linewidth=2)
	plt.text(astronomical_night[-1]-0.7, 0.93, 'Night end', fontsize=10)

	#plot the moon line
	plt.plot(moon_delta_midnight, moon_airmass, 'w-', label='Moon', linestyle='dashed')

	len_range = int(len(nightrange))
	mid_len_range = int(len(nightrange)/2)

	#plot the blue gradients for the astronomical night time
	for xi in range(0, mid_len_range+1):
		plt.fill_between(nightrange[xi:mid_len_range+1], 3, 1, sun_alt_beforemid[xi:mid_len_range+1], facecolor='blue', zorder=0, alpha=0.008)
	for xi in range(len_range, mid_len_range, -1):
		plt.fill_between(nightrange[mid_len_range:xi], 3, 1, sun_alt_beforemid[mid_len_range:xi], facecolor='blue', zorder=0, alpha=0.008)

	#plot the scheduled target info
	for t in scheduled_targets:
		t_altaz = t.skycoord.transform_to(frame_night)
		dm = [x for i,x in enumerate(delta_midnight) if t_altaz.secz[i] > 1 and t_altaz.secz[i] <= 3.5]
		cc = [x for x in t_altaz.secz if x > 1 and x <= 3.5]

		#plot the airmass chart
		plt.plot(dm, cc, label=t.name, linestyle='dashed')

		s_start = (t.start_observation-midnight).to('hour')/u.hour
		s_end = (t.end_observation-midnight).to('hour')/u.hour

		scheduled_time_range = np.linspace(s_start, s_end, 50)
		bool_range = [True for x in scheduled_time_range]
		sc_aa_tmp = [x for i,x in enumerate(cc) if dm[i] > s_start and dm[i] < s_end]
		scheduled_altaz = np.linspace(sc_aa_tmp[0], sc_aa_tmp[-1], 50)

		#plot the observation range as a filled column
		plt.fill_between(scheduled_time_range, 3, 1, bool_range, facecolor='red', zorder=0, alpha=0.4)
		#plot the observation range on top of the airmass line
		plt.plot(scheduled_time_range, scheduled_altaz, color='red', linewidth=10, alpha=0.4)

		plt.text((s_end+s_start)/2.0, 2.625, t.name, ha='center', va='center', rotation=90, size=9)

		print('{}: {} - {}. {} {}'.format(
			t.name,
			round(float(scheduled_time_range[0]), 3),
			round(float(scheduled_time_range[-1]), 3),
			t.priority,
			np.mean(scheduled_altaz)
		))

	#plt.tight_layout()
	fig.subplots_adjust(right=0.75)
	plt.legend(loc='upper left', bbox_to_anchor=(1, 1), borderaxespad=0.)
	plt.grid(True, color='k')
	plt.xlim(nightrange[0], nightrange[len(nightrange)-1])
	plt.ylim(3, 0.95)
	plt.xlabel('Hours from EDT Midnight')
	plt.ylabel('Airmass')
	plt.show()

	#populate_rts2_queue(scheduled_targets)


def populate_rts2_queue(scheduled_targets, set_queue='plan'):
    q = queue.Queue(set_queue)
    print('Clearing current Queue')
    q.load()
    q.clear()
    print('Scheduling')
    for t in scheduled_targets:
        print(t.name)
        rts2ids = t.rts2ids
        for rid in t.rts2ids:
            print('   ', rid)
            q.add_target(rid)

    print('Queue populated')


def oldmain():

	bigelow = astroplan.Observer(latitude=32.4151*astropy.units.deg,
								 longitude=110.7145*astropy.units.deg,
								 elevation=2607*astropy.units.m,
								 name='Bigelow',
								 timezone='America/Phoenix')

	lemmon_sc = astroplan.Observer(latitude=32.4420*astropy.units.deg,
								 longitude=110.7893*astropy.units.deg,
								 elevation=2791*astropy.units.m,
								 name='Mt Lemmon',
								 timezone='America/Phoenix')

	observatories = [bigelow, lemmon_sc]

	targets = readfromweb()

	target_table = [astroplan.FixedTarget(
			coord=astropy.coordinates.SkyCoord(
				ra=astropy.coordinates.Angle(d.ra),
				dec=astropy.coordinates.Angle(d.dec)),
			name=d.name)
	for d in targets]

	for d in targets:
		d.outputobjectinfo()

	global_constraints = [astroplan.constraints.AirmassConstraint(max = 3, boolean_constraint = False),
						  astroplan.constraints.AtNightConstraint.twilight_civil()]
	blocks = []
	read_out = 10*astropy.units.second
	Night = astroplan.constraints.TimeConstraint(astropy.time.Time('2021-03-04 19:00'), astropy.time.Time('2021-03-05 06:00'))

	for t in target_table:
		targetinfo = [x.exp_objs for x in targets if t.name == x.name][0]
		for exp in targetinfo:
			priority = 1
			n = int(exp.amount)
			exptime = int(exp.exptime)*astropy.units.second
			b = astroplan.ObservingBlock.from_exposures(t,
				priority,
				exptime,
				n,
				read_out,
				configuration= {'filter':exp.filter},
				constraints = [Night])
			blocks.append(b)


	#Transitioner
	slew_rate = 0.8*astropy.units.deg/astropy.units.second
	transitioner = astroplan.scheduling.Transitioner(slew_rate,
					{'filter':{'default':10*astropy.units.second}})


	seq_scheduler = astroplan.scheduling.SequentialScheduler(constraints = global_constraints,
										observer = bigelow,
										transitioner = transitioner)

	noon_before = astropy.time.Time('2021-03-04 12:00')
	noon_after = astropy.time.Time('2021-03-05 12:00')

	sequential_schedule = astroplan.scheduling.Schedule(noon_before, noon_after)
	print("this might take a minute")
	seq_scheduler(blocks, sequential_schedule)
	print("yiff")
	print(sequential_schedule.to_table())

main()
