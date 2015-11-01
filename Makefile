all:
	@javac SoundFile.java MatchAudio.java Peak.java
	@chmod 755 dam

clean:
	@rm -rf SoundFile.class MatchAudio.class Peak.class
