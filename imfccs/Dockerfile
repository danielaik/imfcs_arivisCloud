FROM maven:3-jdk-11 as builder

COPY . /usr/src/app
WORKDIR /usr/src/app

RUN mvn clean compile assembly:single

FROM openjdk:11-jre-slim
WORKDIR /usr/src/app

COPY ./template/template_v8_FCCS.html .

COPY --from=builder /usr/src/app/target/FCCS-multitau-correlator-1.0.0-jar-with-dependencies.jar .
CMD [ "java", "-jar", "FCCS-multitau-correlator-1.0.0-jar-with-dependencies.jar" ]
